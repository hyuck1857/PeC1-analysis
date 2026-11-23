"""Check the 31 calibration jobs in Figure S3a against Supplementary Table S5.

Each CSV row is the ipTM of one job's highest-ranked model among five.
Repeated scores are separate jobs. Values come from the iptm column of
C1_interaction_boost_analysis/results/tables/07_calibration_master.csv.
Rows are sorted by score within each group and omit the original job names.

The group counts, arithmetic means and extremes below adapt the study's
C1_interaction_boost_analysis/scripts/06_parse_calibration_results.py,
which wrote 08_calibration_summary.csv. The inclusive counts within the
observed negative-control range are an added descriptive check. Table S5
reports a separate 15-job analysis using a different score summary.
"""

import csv
from decimal import Decimal
from pathlib import Path
from statistics import mean


DATA = Path(__file__).resolve().parents[1] / "data" / "interface_confidence.csv"

# Table S5 reports group size, mean (shown here to three decimals) and range.
TABLE_S5 = {
    "positive_control": (2, "0.180", "0.13", "0.23"),
    "negative_control": (2, "0.120", "0.10", "0.14"),
    "catalytic_domain_panel": (13, "0.168", "0.11", "0.21"),
    "structural_domain_panel": (14, "0.168", "0.08", "0.25"),
}


def main():
    scores_by_group = {group: [] for group in TABLE_S5}
    with DATA.open(encoding="utf-8", newline="") as source:
        for row in csv.DictReader(source):
            group, score = row["group"], Decimal(row["iptm"])
            if group not in TABLE_S5 or not score.is_finite() or not 0 <= score <= 1:
                raise ValueError(f"Invalid calibration row: {row}")
            scores_by_group[group].append(score)

    for group, (expected_n, *_) in TABLE_S5.items():
        if len(scores_by_group[group]) != expected_n:
            raise ValueError(f"Table S5 reports {expected_n} jobs for {group}.")

    negatives = scores_by_group["negative_control"]
    low, high = min(negatives), max(negatives)
    print("group,n,mean_iptm,min_iptm,max_iptm,n_in_negative_range")
    for group, scores in scores_by_group.items():
        summary = (len(scores), f"{mean(scores):.3f}",
                   f"{min(scores):.2f}", f"{max(scores):.2f}")
        if summary != TABLE_S5[group]:
            raise ValueError(f"Summary for {group} differs from Table S5: {summary}")
        overlap = sum(low <= score <= high for score in scores)
        print(group, *summary, overlap, sep=",")


if __name__ == "__main__":
    main()
