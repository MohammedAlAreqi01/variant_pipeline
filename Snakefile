# Snakefile

# gnomAD variant pathogenicity pipeline (CSV-based, reproduces the manual
# thesis workflow end-to-end).
#
# Per gene you drop a raw gnomAD browser CSV export into  input/{gene}.{version}.csv
# and the pipeline runs:
#   raw gnomAD CSV
#     -> avinput            (gnomad_to_avinput.py)
#     -> ANNOVAR annotate   (table_annovar.pl, dbnsfp42a + clinvar)
#     -> InterVar           (kept for future inframe-indel ACMG)
#     -> merge predictors back onto gnomAD columns (merge_annotations.py)
#     -> classify           (classify.py - validated Arm B criteria)
#     -> summary tables      (summarize.py)
#
# Run:  snakemake --cores 4
# One gene:  snakemake --cores 4 results/PALB2.v4.1.summary.txt

configfile: "config.yaml"

GENES = list(config["genes"].keys())
GNOMAD_VERSIONS = config["gnomad_versions"]

wildcard_constraints:
    gene="|".join(GENES),
    version="|".join(v.replace(".", r"\.") for v in GNOMAD_VERSIONS)


def annovar_build(wc):
    return "hg38" if config["gnomad"][wc.version]["build"] == "GRCh38" else "hg19"


rule all:
    input:
        expand("results/{gene}.{version}.classified.csv",
               gene=GENES, version=GNOMAD_VERSIONS),
        expand("results/{gene}.{version}.summary.txt",
               gene=GENES, version=GNOMAD_VERSIONS),
        # InterVar isn't consumed by classify/merge yet (kept for future
        # inframe-indel handling), but is still required here so a plain
        # `snakemake --cores 4` actually runs it instead of silently
        # skipping a rule nothing else depends on.
        expand("intervar/{gene}.{version}.intervar.tsv",
               gene=GENES, version=GNOMAD_VERSIONS)


# ---------------------------------------------------------------------
# 1. gnomAD CSV -> ANNOVAR avinput
# ---------------------------------------------------------------------
rule gnomad_to_avinput:
    input:
        "input/{gene}.{version}.csv"
    output:
        "avinput/{gene}.{version}.avinput"
    log:
        "logs/avinput_{gene}_{version}.log"
    shell:
        """
        mkdir -p avinput
        python scripts/gnomad_to_avinput.py \
            --input {input} --output {output} > {log} 2>&1
        """


# ---------------------------------------------------------------------
# 2. ANNOVAR annotation (refGene, dbnsfp42a, clinvar) - csvout
# ---------------------------------------------------------------------
rule annovar_annotate:
    input:
        "avinput/{gene}.{version}.avinput"
    output:
        "annotated/{gene}.{version}.multianno.csv"
    params:
        annovar_dir = config["tools"]["annovar_dir"],
        humandb = config["tools"]["humandb_dir"],
        build = annovar_build,
        prefix = "annotated/{gene}.{version}"
    log:
        "logs/annovar_{gene}_{version}.log"
    shell:
        """
        mkdir -p annotated
        {params.annovar_dir}/table_annovar.pl {input} {params.humandb} \
            -buildver {params.build} \
            -out {params.prefix} \
            -remove -protocol refGene,dbnsfp42a,clinvar_20250721 \
            -operation g,f,f -nastring . -csvout \
            > {log} 2>&1
        mv {params.prefix}.{params.build}_multianno.csv {output}
        """


# ---------------------------------------------------------------------
# 3. InterVar (ACMG/AMP) - kept for future inframe-indel handling
# ---------------------------------------------------------------------
rule intervar_classify:
    input:
        "avinput/{gene}.{version}.avinput"
    output:
        "intervar/{gene}.{version}.intervar.tsv"
    params:
        intervar_dir = config["tools"]["intervar_dir"],
        humandb = config["tools"]["humandb_dir"],
        build = annovar_build,
        prefix = "intervar/{gene}.{version}"
    log:
        "logs/intervar_{gene}_{version}.log"
    shell:
        """
        mkdir -p intervar
        python {params.intervar_dir}/Intervar.py \
            -i {input} -o {params.prefix} \
            -d {params.humandb} -b {params.build} \
            > {log} 2>&1 || true
        if [ -f {params.prefix}.{params.build}_multianno.txt.intervar ]; then
            mv {params.prefix}.{params.build}_multianno.txt.intervar {output}
        else
            echo "InterVar produced no output" >> {log}
            touch {output}
        fi
        """


# ---------------------------------------------------------------------
# 4. Merge ANNOVAR predictor columns back onto gnomAD columns
# ---------------------------------------------------------------------
rule merge_annotations:
    input:
        gnomad = "input/{gene}.{version}.csv",
        annovar = "annotated/{gene}.{version}.multianno.csv"
    output:
        "merged/{gene}.{version}.merged.csv"
    log:
        "logs/merge_{gene}_{version}.log"
    shell:
        """
        mkdir -p merged
        python scripts/merge_annotations.py \
            --gnomad {input.gnomad} \
            --annovar {input.annovar} \
            --output {output} > {log} 2>&1
        """


# ---------------------------------------------------------------------
# 5. Classify (validated Arm B criteria, per-gene exon exclusion from config)
# ---------------------------------------------------------------------
rule classify_variants:
    input:
        "merged/{gene}.{version}.merged.csv"
    output:
        "results/{gene}.{version}.classified.csv"
    params:
        cutoff_arg = lambda wc: (
            f"--exon-cutoff {config['genes'][wc.gene]['exon_cutoff']}"
            if config["genes"][wc.gene].get("exon_cutoff") is not None else ""
        ),
        t = config["missense_thresholds"],
    log:
        "logs/classify_{gene}_{version}.log"
    shell:
        """
        mkdir -p results
        python scripts/classify.py \
            --input {input} --output {output} \
            --gene {wildcards.gene} {params.cutoff_arg} \
            --max-allele-count {params.t[max_allele_count]} \
            --polyphen-hvar-min {params.t[polyphen_hvar_min]} \
            --sift-max {params.t[sift_max]} \
            --gerp-nr-min {params.t[gerp_nr_min]} \
            > {log} 2>&1
        """


# ---------------------------------------------------------------------
# 6. Summary (Table 1 class counts + Table 2 per-ancestry)
# ---------------------------------------------------------------------
rule summarize:
    input:
        "results/{gene}.{version}.classified.csv"
    output:
        summary = "results/{gene}.{version}.summary.txt",
        ancestry = "results/{gene}.{version}.ancestry.csv"
    log:
        "logs/summarize_{gene}_{version}.log"
    shell:
        """
        python scripts/summarize.py \
            --input {input} \
            --output {output.ancestry} \
            > {output.summary} 2> {log}
        """
