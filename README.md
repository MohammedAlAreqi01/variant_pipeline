# gnomAD Variant Pathogenicity Pipeline

An automated Snakemake pipeline that classifies gnomAD variants as
"predicted pathogenic" and reports their population frequencies per genetic
ancestry group. It reproduces a manual ANNOVAR/InterVar/ACMG workflow,
end-to-end, from a single raw gnomAD CSV export.

Validated against manual PALB2 analysis: reproduces the manual calls exactly
(380 predicted pathogenic; null 220, missense 160), including per-ancestry
cohort sizes.

## What it does

Starting from a raw gnomAD browser CSV export for one gene, the pipeline:

1. Converts the CSV to ANNOVAR input format (avinput)
2. Annotates with ANNOVAR (refGene, dbNSFP 4.2a predictors, ClinVar)
3. Runs InterVar (ACMG/AMP; retained for future inframe-indel handling)
4. Merges ANNOVAR's predictor columns back onto the gnomAD columns
5. Classifies each variant as predicted pathogenic (criteria below)
6. Summarises counts by variant class and by ancestry

## Classification criteria

A variant is "predicted pathogenic" (Arm B, computational prediction) if,
after removing synonymous and non-coding variants, it meets its class rule:

- **Null** (frameshift, stop_gained, splice acceptor/donor): predicted
  pathogenic regardless of allele count, EXCEPT variants in the last two
  exons (NMD escape), which are excluded via a per-gene position cutoff.
- **Missense**: Allele Count <= 5 AND Polyphen2_HVAR > 0.95 AND SIFT < 0.05
  AND MutationTaster prediction in {A, D} AND GERP++_NR >= 3.
- **Inframe indels**: ACMG criteria (not yet automated).
- **Structural / CNV**: separate gnomAD SV dataset (not yet integrated).

"Number of people" for a class = sum of that class's allele counts.
Per-ancestry cohort size = mean Allele Number across all variants in the
ancestry, divided by 2.

## Folder layout

```
Snakefile               the workflow (7 rules)
config.yaml             gene list, gnomAD versions, tool paths
scripts/
  gnomad_to_avinput.py    gnomAD CSV -> ANNOVAR avinput
  merge_annotations.py     join ANNOVAR predictors onto gnomAD columns
  classify.py              apply the classification criteria
  summarize.py             class + per-ancestry summary tables
input/                  drop raw gnomAD CSVs here as {gene}.{version}.csv
```

Generated at run time (safe to delete; regenerated): `avinput/`,
`annotated/`, `intervar/`, `merged/`, `results/`, `logs/`.

## Requirements

- Conda (Miniconda/Mambaforge). Create the environment from the included file:
  ```bash
  conda env create -f environment.yml
  conda activate snakemake
  ```
  This installs Snakemake, tabix/bgzip, and bcftools. The pipeline's own
  scripts use only the Python standard library, so nothing else is needed.
- ANNOVAR with the `humandb` databases: refGene, dbnsfp42a, clinvar_20250721
  (register and download separately; not conda-installable)
- InterVar 2.2.1 (git clone; edit its `config.ini` to point at the ANNOVAR install)
- Set the three tool paths in `config.yaml` under `tools:` to match this machine

## Running it

Place a raw gnomAD CSV export at `input/{gene}.{version}.csv`, e.g.
`input/PALB2.v4.1.csv`, then:

```bash
conda activate snakemake

snakemake --cores 4 -n                              # dry run (check the plan)
snakemake --cores 4 results/PALB2.v4.1.summary.txt  # run one gene
snakemake --cores 4                                 # run everything in config
```

Results land in `results/`:
- `{gene}.{version}.classified.csv` - every variant + a predicted_pathogenic
  flag and the rule (arm) that classified it
- `{gene}.{version}.summary.txt` - overall n, class counts, per-ancestry table
- `{gene}.{version}.ancestry.csv` - the per-ancestry table as CSV

## Adding a new gene

1. Export the gene's variants from the gnomAD browser as CSV.
2. Save it as `input/{GENE}.{version}.csv`.
3. Add the gene to `config.yaml` under `genes:`.
4. Add the gene's last-two-exon cutoff to `EXON_EXCLUSION` in
   `scripts/classify.py` (strand-aware GRCh38 position boundary), or its
   null variants in the final exons won't be excluded.

## Status / not yet automated

- Inframe-indel ACMG classification (none qualified for PALB2, but needed
  before extending to genes where they might)
- Structural / CNV variants (separate gnomAD SV dataset, with cohort-size
  correction factors)
- Formatted thesis-style Table 1 / Table 2 output
