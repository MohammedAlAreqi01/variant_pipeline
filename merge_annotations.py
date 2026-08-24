#!/usr/bin/env python3
"""
Merge ANNOVAR's -csvout predictor columns back onto the original gnomAD CSV,
reproducing the (previously manual, in-Excel) combined file that classify.py
consumes.

Join key: chromosome / position / ref / alt.
  gnomAD columns : Chromosome, Position, Reference, Alternate
  ANNOVAR columns: Chr, Start, Ref, Alt

Output = every gnomAD column, followed by the ANNOVAR predictor columns
(the dbNSFP block: SIFT, Polyphen2, MutationTaster, GERP++, etc.), matching
the validated PALB2_v4.1.1.csv layout (gnomAD cols 0-73, ANNOVAR cols 74+).

ANNOVAR structural columns (Chr/Start/End/Ref/Alt/Func.refGene/Gene.refGene)
are used for the join and dropped from the appended block, since the gnomAD
side already carries position info and classify.py keys off the gnomAD/dbNSFP
column names.
"""
import argparse
import csv

# ANNOVAR columns to use for the join but NOT append (avoid duplicating
# position/gene info the gnomAD side already has). Everything else ANNOVAR
# adds (the dbNSFP predictor block) is appended.
ANNOVAR_KEY_COLS = {"Chr", "Start", "End", "Ref", "Alt"}


def norm(chrom):
    """Normalise chromosome to bare form (strip any 'chr' prefix)."""
    return str(chrom).strip().replace("chr", "").replace("Chr", "")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gnomad", required=True, help="original gnomAD CSV")
    p.add_argument("--annovar", required=True, help="ANNOVAR .hg38_multianno.csv")
    p.add_argument("--output", required=True, help="merged CSV for classify.py")
    args = p.parse_args()

    # Read gnomAD CSV
    with open(args.gnomad, newline="", encoding="utf-8-sig") as f:
        greader = csv.reader(f)
        gheader = next(greader)
        grows = list(greader)
    gi = {name: idx for idx, name in enumerate(gheader)}

    # Read ANNOVAR csvout
    with open(args.annovar, newline="", encoding="utf-8-sig") as f:
        areader = csv.reader(f)
        aheader = next(areader)
        arows = list(areader)
    ai = {name: idx for idx, name in enumerate(aheader)}

    # Which ANNOVAR columns to append (predictor block = all except keys)
    append_cols = [name for name in aheader if name not in ANNOVAR_KEY_COLS]
    append_idx = [ai[name] for name in append_cols]

    # Build an ANNOVAR lookup keyed by (chrom, start, ref, alt)
    a_lookup = {}
    for row in arows:
        key = (
            norm(row[ai["Chr"]]),
            str(row[ai["Start"]]).strip(),
            str(row[ai["Ref"]]).strip(),
            str(row[ai["Alt"]]).strip(),
        )
        a_lookup[key] = row

    # gnomAD key columns
    g_chr = gi["Chromosome"]
    g_pos = gi["Position"]
    g_ref = gi["Reference"]
    g_alt = gi["Alternate"]

    out_header = gheader + append_cols
    n_matched = 0
    n_missed = 0
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(out_header)
        for grow in grows:
            key = (
                norm(grow[g_chr]),
                str(grow[g_pos]).strip(),
                str(grow[g_ref]).strip(),
                str(grow[g_alt]).strip(),
            )
            arow = a_lookup.get(key)
            if arow is not None:
                n_matched += 1
                appended = [arow[i] for i in append_idx]
            else:
                n_missed += 1
                appended = ["" for _ in append_idx]  # no annotation match
            writer.writerow(grow + appended)

    print(f"Merged {len(grows)} gnomAD variants: "
          f"{n_matched} matched ANNOVAR, {n_missed} unmatched")


if __name__ == "__main__":
    main()
