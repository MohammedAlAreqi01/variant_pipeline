#!/usr/bin/env python3
"""
Summarise a classified gnomAD CSV (output of classify.py).

Produces:
  1. Breakdown of predicted-pathogenic variants by class (null / missense).
  2. Per-ancestry tallies for predicted-pathogenic variants:
       - variant count  : number of distinct predicted-pathogenic variants
                           observed (Allele Count > 0) in that ancestry
       - allele count    : total predicted-pathogenic alleles in that ancestry
       - allele number   : total allele number (for frequency denominators)
       - frequency       : alleles / allele number
"""
import argparse
import csv

ANCESTRIES = [
    "African/African American", "Admixed American", "Ashkenazi Jewish",
    "East Asian", "European (Finnish)", "Middle Eastern",
    "European (non-Finnish)", "Amish", "South Asian", "Remaining",
]


def to_int(val):
    val = (val or "").strip()
    if val in ("", ".", "NA"):
        return 0
    try:
        return int(float(val))
    except ValueError:
        return 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="classified CSV from classify.py")
    p.add_argument("--output", help="optional: write per-ancestry table to this CSV")
    args = p.parse_args()

    with open(args.input, newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)

    def idx(name):
        return header.index(name)

    pp_i = idx("predicted_pathogenic")
    arm_i = idx("arm")

    # keep only predicted-pathogenic rows
    pp_rows = [r for r in rows if len(r) > pp_i and r[pp_i] == "yes"]

    # 1. class breakdown (null vs missense, from the 'arm' reason)
    #    with people = sum of total Allele Count within the class
    from collections import Counter
    ac_total_i = idx("Allele Count")
    by_arm_variants = Counter()
    by_arm_people = Counter()
    for r in pp_rows:
        by_arm_variants[r[arm_i]] += 1
        by_arm_people[r[arm_i]] += to_int(r[ac_total_i])

    total_people = sum(by_arm_people.values())

    # overall cohort n: mean of the total 'Allele Number' column across all
    # variants (all ancestries combined), divided by 2
    total_an_i = idx("Allele Number")
    total_an_values = [to_int(r[total_an_i]) for r in rows if len(r) > total_an_i]
    total_an_values = [v for v in total_an_values if v > 0]
    overall_n = (sum(total_an_values) / len(total_an_values) / 2) if total_an_values else 0

    print("=== Predicted-pathogenic by class ===")
    print(f"  Overall cohort n (avg Allele Number / 2): {overall_n:,.0f}")
    print(f"  Total predicted pathogenic : {len(pp_rows)} variants / {total_people} people")
    for arm, n in by_arm_variants.most_common():
        label = arm.replace("ArmB:", "") if arm else "(unlabelled)"
        print(f"  {label:12s}: {n} variants / {by_arm_people[arm]} people")
    print()

    # 2. per-ancestry tallies
    #    variant count : distinct PP variants with AC>0 in the ancestry
    #    people        : sum of PP allele counts in the ancestry
    #    cohort size   : mean Allele Number across ALL variants in the
    #                    ancestry, divided by 2 (AN is ~2x people)
    #    one-in-X       : cohort size / people
    print("=== Predicted-pathogenic per ancestry ===")
    print(f"  {'Ancestry':26s} {'Variants':>9s} {'People':>8s} {'Cohort':>10s} {'1 in':>9s}")
    out_rows = []
    for anc in ANCESTRIES:
        ac_i = idx(f"Allele Count {anc}")
        an_i = idx(f"Allele Number {anc}")

        # cohort size: average AN across ALL rows (not just PP), halved
        an_values = [to_int(r[an_i]) for r in rows if len(r) > an_i]
        an_values = [v for v in an_values if v > 0]
        mean_an = (sum(an_values) / len(an_values)) if an_values else 0
        cohort = mean_an / 2

        # variant count + people on the PP set
        variant_count = 0
        people = 0
        for r in pp_rows:
            ac = to_int(r[ac_i])
            if ac > 0:
                variant_count += 1
                people += ac

        one_in = (cohort / people) if people else 0
        one_in_str = f"1 in {one_in:,.0f}" if people else "-"
        print(f"  {anc:26s} {variant_count:9d} {people:8d} {cohort:10,.0f} {one_in_str:>12s}")
        out_rows.append({
            "ancestry": anc,
            "predicted_pathogenic_variants": variant_count,
            "people_alleles": people,
            "cohort_size": round(cohort),
            "one_in": round(one_in) if people else "",
        })

    if args.output:
        with open(args.output, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
            w.writeheader()
            w.writerows(out_rows)
        print(f"\nPer-ancestry table written to {args.output}")


if __name__ == "__main__":
    main()
