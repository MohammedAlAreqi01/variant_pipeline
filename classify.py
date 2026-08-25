#!/usr/bin/env python3
"""
Classify gnomAD variants as "Predicted pathogenic" per Mo's thesis criteria.

Two arms (see thesis flowchart):

  ARM A - ClinVar:
    Variants classified Pathogenic or Likely pathogenic in the gnomAD-provided
    'ClinVar Germline Classification' column.

  ARM B - Computationally predicted (after removing synonymous & non-coding):
    - Missense: Allele Count <= 5 AND Polyphen2_HVAR_score > 0.95
      AND SIFT_score < 0.05 AND MutationTaster_pred in {A, D} AND GERP++_NR >= 3
    - Null (frameshift, stop_gained, splice acceptor/donor): Allele Count <= 5
      [Pass 2 TODO: exclude last 50bp / penultimate exon]
    - Structural (DEL/DUP): LoF only, Allele Count <= 5
      [not represented as SVs in the gnomAD short-variant CSV export]
    - Inframe indels: ACMG criteria (via InterVar)  [Pass 2 TODO]

A variant is "Predicted pathogenic" if it qualifies under EITHER arm.

Input: the per-gene gnomAD CSV export (with ANNOVAR predictor columns appended).
Output: same rows plus a 'predicted_pathogenic' flag and 'arm' reason column.
"""
import argparse
import csv
import re
import sys


# VEP consequence -> variant class (col 14 "VEP Annotation")
# Null / LoF variants. Per validation against the thesis calls, these are
# classified predicted-pathogenic regardless of allele count (the AC<=5
# filter applies ONLY to missense). start_lost and protein_altering_variant
# are deliberately NOT included.
NULL_CONSEQUENCES = {
    "frameshift_variant", "stop_gained",
    "splice_acceptor_variant", "splice_donor_variant",
}
MISSENSE_CONSEQUENCES = {"missense_variant"}
INFRAME_CONSEQUENCES = {"inframe_deletion", "inframe_insertion"}
# Explicitly excluded (synonymous & non-coding)
NONCODING_CONSEQUENCES = {
    "intron_variant", "synonymous_variant",
    "5_prime_UTR_variant", "3_prime_UTR_variant",
    "splice_region_variant",  # non-LoF splice region; not a null variant
}

CLINVAR_PATHOGENIC = {
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
}


# Pass 2: penultimate-exon exclusion for null variants.
# Null variants beyond the last 50bp of the penultimate exon escape
# nonsense-mediated decay and are excluded (per thesis criteria).
#
# The boundary is a per-gene cDNA (HGVS c.) position from Alamut - the last
# 50bp of the penultimate exon - e.g. PALB2 c.3300, ENG c.1802, ACVRL1 c.1327
# (MANE Select transcripts). A null variant is excluded when its coding
# position is GREATER THAN the cutoff (i.e. 3' of it, toward the C-terminus).
#
# Filtering on c. coordinates (not genomic) matches the manual Alamut method
# and works identically regardless of gene strand, since c. numbering already
# runs 5'->3'.

def parse_c_pos(hgvs):
    """Extract the coding (c.) position from an HGVS Transcript Consequence.

    Returns the highest coding position in the variant (closest to the 3' end,
    so range variants spanning the boundary are judged by their 3'-most base),
    ignoring intronic +N/-N offsets. Returns None if unparseable.

    Examples:
      c.3340C>T          -> 3340
      c.3295_3305del     -> 3305
      c.3350+14A>G       -> 3350   (intronic offset stripped)
      c.3256_3257ins...  -> 3257
    """
    if not hgvs or not str(hgvs).startswith("c."):
        return None
    s = str(hgvs)[2:]
    # coding positions are integers NOT immediately preceded by + or - (those
    # mark intronic offsets like 3350+14, whose base coding position is 3350)
    positions = [int(m.group(1))
                 for m in re.finditer(r"(?<![+\-\d])(\d+)", s)]
    return max(positions) if positions else None


def in_excluded_exon(hgvs_c, cutoff):
    """True if a null variant lies 3' of the penultimate-exon cutoff.

    cutoff is the gene's c. boundary (e.g. 3300 for PALB2). A variant is
    excluded when its coding position is strictly greater than the cutoff.
    If no cutoff is configured (None), nothing is excluded.
    """
    if cutoff is None:
        return False
    c_pos = parse_c_pos(hgvs_c)
    if c_pos is None:
        return False
    return c_pos > cutoff



def to_float(val):
    """Parse a numeric cell; return None if blank or non-numeric ('.' etc.)."""
    if val is None:
        return None
    val = val.strip()
    if val in ("", ".", "NA", "na"):
        return None
    # some dbNSFP cells pack multiple transcript scores like "0.99;0.87"
    if ";" in val:
        parts = [p for p in val.split(";") if p not in ("", ".")]
        if not parts:
            return None
        try:
            # use the max (most damaging for score-higher-is-worse) - but
            # for SIFT lower is worse, so caller decides; return max here and
            # handle SIFT specially below
            return max(float(p) for p in parts)
        except ValueError:
            return None
    try:
        return float(val)
    except ValueError:
        return None


def to_float_min(val):
    """Like to_float but returns the MIN across multi-transcript cells (for SIFT)."""
    if val is None:
        return None
    val = val.strip()
    if val in ("", ".", "NA", "na"):
        return None
    if ";" in val:
        parts = [p for p in val.split(";") if p not in ("", ".")]
        if not parts:
            return None
        try:
            return min(float(p) for p in parts)
        except ValueError:
            return None
    try:
        return float(val)
    except ValueError:
        return None


def mt_is_disease_causing(pred):
    """MutationTaster prediction is 'A' or 'D' (possibly multi-transcript)."""
    if not pred:
        return False
    tokens = pred.replace(";", " ").split()
    return any(t in ("A", "D") for t in tokens)


def classify_row(row, cols):
    """Return (is_pathogenic, arm_reason) for one variant row.

    Arm B only (computational prediction). The ClinVar arm is intentionally
    not applied - classification is purely the consequence + AC + predictor
    filters below.
    """
    vep = row[cols["vep"]].strip()

    # --- ARM B: remove synonymous & non-coding up front ---
    if vep in NONCODING_CONSEQUENCES or vep == "":
        return False, ""

    # Null variants: predicted-pathogenic regardless of allele count,
    # EXCEPT those past the penultimate-exon c. cutoff (NMD escape).
    if vep in NULL_CONSEQUENCES:
        if in_excluded_exon(row[cols["hgvs_c"]], cols["exon_cutoff"]):
            return False, "ArmB:null_excluded_exon"
        return True, "ArmB:null"

    # Missense rule: AC<=cutoff AND all four predictor criteria.
    # Thresholds come from config (with defaults matching the thesis).
    if vep in MISSENSE_CONSEQUENCES:
        t = cols["thresholds"]
        ac = to_float(row[cols["ac"]])
        if ac is None or ac > t["max_allele_count"]:
            return False, ""
        polyphen = to_float(row[cols["polyphen_hvar"]])
        sift = to_float_min(row[cols["sift"]])
        mt = row[cols["mt_pred"]]
        gerp = to_float(row[cols["gerp_nr"]])
        if (polyphen is not None and polyphen > t["polyphen_hvar_min"]
                and sift is not None and sift < t["sift_max"]
                and mt_is_disease_causing(mt)
                and gerp is not None and gerp >= t["gerp_nr_min"]):
            return True, "ArmB:missense"
        return False, ""

    # Inframe indels - Pass 2 TODO: ACMG criteria via InterVar.
    if vep in INFRAME_CONSEQUENCES:
        return False, "inframe_pending_ACMG"

    return False, ""


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--gene", required=True,
                   help="gene symbol (e.g. PALB2)")
    p.add_argument("--exon-cutoff", type=int, default=None,
                   help="cDNA (c.) position cutoff for penultimate-exon "
                        "exclusion, from Alamut (e.g. 3300 for PALB2). Null "
                        "variants with c.position > cutoff are excluded. "
                        "Omit to disable exclusion.")
    # Missense thresholds (defaults match the validated thesis methodology).
    # Supplied from config.yaml via the Snakefile so they can be tuned or
    # extended (e.g. adding AlphaMissense) without editing this script.
    p.add_argument("--max-allele-count", type=float, default=5,
                   help="missense: max Allele Count (default 5)")
    p.add_argument("--polyphen-hvar-min", type=float, default=0.95,
                   help="missense: Polyphen2_HVAR must be > this (default 0.95)")
    p.add_argument("--sift-max", type=float, default=0.05,
                   help="missense: SIFT must be < this (default 0.05)")
    p.add_argument("--gerp-nr-min", type=float, default=3,
                   help="missense: GERP++_NR must be >= this (default 3)")
    args = p.parse_args()

    with open(args.input, newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)

    def find(name):
        try:
            return header.index(name)
        except ValueError:
            sys.exit(f"ERROR: expected column {name!r} not found in input")

    def find_any(*names):
        """Return the index of the first column name that exists.

        Accepts either ANNOVAR's real output names (e.g. 'SIFT_score') or the
        threshold-encoded names from the original hand-edited thesis file
        (e.g. 'SIFT_score<0.05'), so the classifier works on both.
        """
        for n in names:
            if n in header:
                return header.index(n)
        sys.exit(f"ERROR: none of these columns found: {names}")

    cols = {
        "vep": find("VEP Annotation"),
        # Looked up but not read in classify_row: Arm A (ClinVar) is
        # intentionally not applied yet (see module docstring). Kept here
        # so it's a one-line change to wire back in.
        "clinvar": find("ClinVar Germline Classification"),
        "ac": find("Allele Count"),
        "polyphen_hvar": find_any("Polyphen2_HVAR_score", "Polyphen2_HVAR_score>0.95"),
        "sift": find_any("SIFT_score", "SIFT_score<0.05"),
        "mt_pred": find("MutationTaster_pred"),
        "gerp_nr": find("GERP++_NR"),
        "pos": find("Position"),
        "hgvs_c": find_any("Transcript Consequence", "HGVS Consequence"),
        "gene": args.gene,
        "exon_cutoff": args.exon_cutoff,
        "thresholds": {
            "max_allele_count": args.max_allele_count,
            "polyphen_hvar_min": args.polyphen_hvar_min,
            "sift_max": args.sift_max,
            "gerp_nr_min": args.gerp_nr_min,
        },
    }

    out_header = header + ["predicted_pathogenic", "arm"]
    n_path = 0
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(out_header)
        for row in rows:
            # pad short rows defensively
            if len(row) < len(header):
                row = row + [""] * (len(header) - len(row))
            is_path, arm = classify_row(row, cols)
            if is_path:
                n_path += 1
            writer.writerow(row + ["yes" if is_path else "no", arm])

    print(f"Classified {len(rows)} variants; {n_path} predicted pathogenic")


if __name__ == "__main__":
    main()
