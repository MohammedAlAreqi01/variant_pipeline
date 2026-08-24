#!/usr/bin/env python3
"""
Ancestry-stratified comparisons on classified variants: Fisher's exact test
where cell counts are small, Chi-squared with Yates correction otherwise.

Called by the `ancestry_stats` Snakemake rule.
"""
import argparse
import pandas as pd
from scipy.stats import fisher_exact, chi2_contingency


def compare_ancestries(df, ancestry_col="ancestry", group_col="classification"):
    results = []
    table = pd.crosstab(df[ancestry_col], df[group_col])

    # Small-count cells -> Fisher's exact; otherwise Chi-squared w/ Yates
    if (table.values < 5).any():
        # Fisher's exact only handles 2x2 natively - collapse if needed
        stat, pval = fisher_exact(table.values[:2, :2])
        test_used = "fisher_exact"
    else:
        stat, pval, _, _ = chi2_contingency(table.values, correction=True)
        test_used = "chi2_yates"

    results.append({"test": test_used, "statistic": stat, "p_value": pval})
    return pd.DataFrame(results)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    result = compare_ancestries(df)
    result.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
