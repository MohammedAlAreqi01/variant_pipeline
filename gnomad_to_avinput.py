#!/usr/bin/env python3
"""
Convert a gnomAD browser CSV export into an ANNOVAR avinput file.

Faithful to Mo's original convert_gnomad_to_avinput.py - the only change is
that input/output paths are CLI arguments (so it runs per-gene in the
pipeline) instead of being hardcoded.

avinput format (tab-separated): chrom  start  end  ref  alt
  - SNV (len(ref)==len(alt)) : end = pos
  - deletion (len(ref)>len(alt)) : end = pos + len(ref) - 2
  - insertion (else)             : end = pos
"""
import argparse
import csv


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="raw gnomAD CSV export")
    p.add_argument("--output", required=True, help="avinput path to write")
    args = p.parse_args()

    with open(args.input, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    with open(args.output, "w", encoding="utf-8") as out:
        for row in rows:
            chrom = str(row["Chromosome"]).strip()
            pos = int(row["Position"])
            ref = str(row["Reference"]).strip()
            alt = str(row["Alternate"]).strip()
            if len(ref) == len(alt):
                end = pos
            elif len(ref) > len(alt):
                end = pos + len(ref) - 2
            else:
                end = pos
            out.write(f"{chrom}\t{pos}\t{end}\t{ref}\t{alt}\n")

    print(f"Done! Wrote avinput to {args.output} ({len(rows)} variants)")


if __name__ == "__main__":
    main()
