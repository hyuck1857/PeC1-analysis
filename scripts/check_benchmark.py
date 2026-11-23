"""Retrospectively check reported phage-level calls against published references.

The accompanying CSV transcribes manuscript Table 2. The five sensitivity
calls come from Supplementary Data S1, sheet 'Class scores', rows marked
'Selected for phage' in the 'T3 sensitivity' condition. The original and
control calls were also cross-checked against the selected rows in that sheet.
Blank accepted_classes means ungraded. The pipe separates accepted classes.
Class labels follow the workbook, where 'colanic' means colanic acid.

This script checks the scorecard from existing calls and performs no prediction.
Run with Python 3 from any directory. No external packages are required.
"""

import csv
from fractions import Fraction
from pathlib import Path


CALLS = Path(__file__).resolve().parents[1] / "data" / "benchmark_calls.csv"
CLASSES = ("cellulose", "LPS", "colanic", "EPS")
TABLE_2_TOTALS = {
    "tier1_control": 3,  # Answer supplied; process control only.
    "tier2": 3,
    "tier3_original": 2,
    "tier3_control": 2,
    "tier3_sensitivity": 1,
}


def main():
    with CALLS.open(encoding="utf-8", newline="") as source:
        rows = list(csv.DictReader(source))
    if len(rows) != 5 or len({row["phage"] for row in rows}) != len(rows):
        raise ValueError("Expected five distinct study phages.")

    graded = [(row, set(row["accepted_classes"].split("|")))
              for row in rows if row["accepted_classes"]]
    if len(graded) != 3:
        raise ValueError("Expected three phages with experimental reference classes.")
    for row, accepted in graded:
        if not accepted.issubset(CLASSES):
            raise ValueError(f"Unknown reference class for {row['phage']}.")
        if any(row[condition] not in CLASSES for condition in TABLE_2_TOTALS):
            raise ValueError(f"Missing or unknown call for {row['phage']}.")

    print(f"Study phages: {len(rows)}; graded: {len(graded)}")
    for condition, expected in TABLE_2_TOTALS.items():
        correct = sum(row[condition] in accepted for row, accepted in graded)
        print(f"{condition}: {correct}/{len(graded)}")
        if correct != expected:
            raise ValueError(f"{condition} differs from Table 2 ({expected} correct).")

    fixed = {label: sum(label in accepted for _, accepted in graded)
             for label in CLASSES}
    best = max(fixed.values())
    best_labels = [label for label, count in fixed.items() if count == best]
    random_expectation = sum(Fraction(len(accepted), len(CLASSES))
                             for _, accepted in graded)
    print(f"Best fixed class ({', '.join(best_labels)}): {best}/{len(graded)}")
    print(f"Expected correct calls under uniform random assignment: {random_expectation}")
    if best != 2 or random_expectation != 1:
        raise ValueError("Calculated baselines differ from Table 2.")


if __name__ == "__main__":
    main()
