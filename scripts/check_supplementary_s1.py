"""Audit the numerical tables in the deposited, final Supplementary Data S1.

Requires Python 3.9+, openpyxl and NumPy. Run from any directory. All comparisons
use archived results from this workbook and the published tables. Counts
printed at the end refer to checks that passed. No external input is needed.
"""

import csv
import hashlib
from collections import Counter, defaultdict
from math import isclose
from pathlib import Path
from statistics import mean, pstdev

import numpy as np
from openpyxl import load_workbook


DATA = Path(__file__).resolve().parents[1] / "data"
WORKBOOK = DATA / "Supplementary_data_S1_v16_20261102.xlsx"
SOURCE_SHA256 = "77bad6e34c399582379f1b5ecd62bc5d480f6f41d616b19735d58df0d381f64e"

# The three columns describe the same conditions in the two S1 sheets
# and the repository's Table 2 scorecard, respectively.
CONDITIONS = {
    "T3 original": ("Original class", "Original 10 class", "tier3_original"),
    "T3 control": ("Control class", "Control 14 class", "tier3_control"),
    "T3 sensitivity": ("Sensitivity class", "Sensitivity 14 class", "tier3_sensitivity"),
}
CLASS_LABELS = {"cellulose": "Cellulose", "LPS": "LPS", "colanic": "Colanic acid", "EPS": "EPS"}

# Table S8, in print order. Each cell triple is (called class, margin, selected).
# The two triples represent the original and expanded control panels.
TABLE_S8 = {
    "phiEaP-21_0042": (("EPS", 0.16, "No"), ("EPS", 0.81, "Yes")),
    "phiEaP-21_0029": (("EPS", 0.46, "No"), ("LPS", 0.05, "No")),
    "phiEaP-21_0036": (("cellulose", 0.56, "Yes"), ("cellulose", 0.49, "No")),
    "phiEaP-21_0046": (("cellulose", 0.63, "No"), ("cellulose", 0.54, "No")),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check_number(actual, recorded, location):
    if not isclose(actual, recorded, abs_tol=1e-10, rel_tol=0):
        raise ValueError(f"{location}: calculated {actual:.15g}, recorded {recorded:.15g}")


def read_sheet(workbook, sheet_name, count):
    # All four S1 tables have column headings on row 5 and data from row 6.
    records = workbook[sheet_name].iter_rows(min_row=5, values_only=True)
    headers = next(records)
    rows = [dict(zip(headers, row)) for row in records]
    require(len(rows) == count, f"{sheet_name}: expected {count} data rows")
    return rows


def read_csv(name):
    with (DATA / name).open(encoding="utf-8", newline="") as source:
        return list(csv.DictReader(source))


def check_ligand_tables(runs, scores, assignments):
    by_ligand = {row["Ligand"]: row for row in assignments}
    require(len(by_ligand) == 14, "Ligand classes: expected 14 distinct ligands")
    batches = Counter(row["Batch"] for row in runs)
    require(batches == {"Original 10": 880, "Added 4": 352}, "Raw runs: unexpected batch sizes")

    seeds = defaultdict(dict)
    for row in runs:
        pair = (row["Candidate RBP"], row["Ligand"])
        seed = row["Seed"]
        require(seed not in seeds[pair], f"Raw runs: duplicate seed {seed} for {pair}")
        seeds[pair][seed] = row["Docking score (kcal mol-1)"]
        require(row["Batch"] == by_ligand[pair[1]]["Panel membership"], f"Raw runs: {pair} batch")

    scores_by_ligand = defaultdict(list)
    for row in scores:
        pair = (row["Candidate RBP"], row["Ligand"])
        values = seeds[pair]
        require(set(values) == set(range(1, 9)), f"Raw runs: incomplete seeds for {pair}")
        best = min(values.values())
        check_number(best, row["Best docking score (kcal mol-1)"], f"Ligand scores: {pair} minimum")
        scores_by_ligand[pair[1]].append(best)
        mapping = by_ligand[pair[1]]
        require(row["Panel membership"] == mapping["Panel membership"], f"Ligand scores: {pair} panel")
        for score_column, mapping_column, _ in CONDITIONS.values():
            require(row[score_column] == mapping[mapping_column], f"Ligand classes: {pair} {mapping_column}")

    z_scores = {}
    for row in scores:
        pair = (row["Candidate RBP"], row["Ligand"])
        sample = scores_by_ligand[pair[1]]
        require(len(sample) == 11 and pair not in z_scores, f"Ligand scores: invalid group for {pair}")
        # The workbook standardizes the negatives of the archived minima.
        z = (mean(sample) - min(seeds[pair].values())) / pstdev(sample)
        check_number(z, row["Docking z-score"], f"Ligand scores: {pair} z-score")
        z_scores[pair] = z
    return z_scores, by_ligand


def check_class_tables(rows, z_scores, assignments):
    groups = defaultdict(list)
    by_candidate = {}
    for row in rows:
        condition, candidate = row["Condition"], row["Candidate RBP"]
        location = f"Class scores: {condition}, {candidate}"
        _, mapping_column, _ = CONDITIONS[condition]
        for label in CLASS_LABELS.values():
            z = [z_scores[candidate, ligand] for ligand, fields in assignments.items()
                 if fields[mapping_column] == label]
            check_number(mean(z), row[label + " score"], f"{location}, {label}")

        called = CLASS_LABELS[row["Predicted class"]]
        other_scores = [row[label + " score"] for label in CLASS_LABELS.values() if label != called]
        difference = row[called + " score"] - max(other_scores)
        require(difference > 0, f"{location}: recorded class does not have a unique highest score")
        check_number(difference, row["Margin"], f"{location}, margin")
        groups[condition, row["Phage"]].append(row)
        key = (condition, candidate)
        require(key not in by_candidate, f"{location}: duplicate row")
        by_candidate[key] = row

    benchmark = read_csv("benchmark_calls.csv")
    expected_groups = {(condition, b["phage"]) for condition in CONDITIONS for b in benchmark}
    require(set(groups) == expected_groups and len(groups) == 15, "Class scores: phage/condition coverage")
    for condition, (_, _, benchmark_column) in CONDITIONS.items():
        for b in benchmark:
            group = groups[condition, b["phage"]]
            selected = [row for row in group if row["Selected for phage"] == "Yes"]
            require(len(selected) == 1, f"{condition}, {b['phage']}: expected one selection")
            chosen = selected[0]
            # All five archived groups contain at least one candidate at pLDDT >= 50.
            eligible = [row for row in group if row["pLDDT"] >= 50]
            require(chosen in eligible, f"{condition}, {b['phage']}: selected pLDDT below 50")
            require(all(chosen["Margin"] > row["Margin"] for row in eligible if row is not chosen),
                    f"{condition}, {b['phage']}: selected margin is not highest among eligible rows")
            require(chosen["Predicted class"] == b[benchmark_column],
                    f"{condition}, {b['phage']}: selected class differs from the benchmark scorecard")
    return by_candidate


def check_supplementary_tables(class_rows, z_scores):
    plddt, matched_means = [], []
    s7 = read_csv("model_confidence.csv")
    require(len({row["candidate"] for row in s7}) == len(s7) == 11, "Table S7: candidate coverage")
    for row in s7:
        candidate = row["candidate"]
        glucan = z_scores[candidate, "cellotriose"]
        galactan = z_scores[candidate, "Gal-b13-Gal-b13-Gal"]
        for value, column in ((glucan, "z_cellotriose"), (galactan, "z_galactan"),
                              (galactan - glucan, "reported_delta_z")):
            check_number(round(value, 2), float(row[column]), f"Table S7: {candidate} {column}")
        confidence = float(row["local_plddt"])
        for condition in CONDITIONS:
            check_number(class_rows[condition, candidate]["pLDDT"], confidence, f"Table S7: {candidate} pLDDT")
        plddt.append(confidence)
        matched_means.append(mean((glucan, galactan)))
    coefficient = float(np.corrcoef(plddt, matched_means)[0, 1])
    check_number(round(coefficient, 2), -0.54, "Figure S3b: Pearson r")

    for candidate, results in TABLE_S8.items():
        for condition, expected in zip(("T3 original", "T3 control"), results):
            row = class_rows[condition, candidate]
            reported = (row["Predicted class"], round(row["Margin"], 2), row["Selected for phage"])
            require(reported == expected, f"Table S8: {candidate}, {condition}: {reported} != {expected}")
    return coefficient


def main():
    require(hashlib.sha256(WORKBOOK.read_bytes()).hexdigest() == SOURCE_SHA256,
            "Supplementary Data S1 differs from the deposited final workbook")
    workbook = load_workbook(WORKBOOK, read_only=True, data_only=True)
    try:
        runs = read_sheet(workbook, "Raw runs", 1232)
        scores = read_sheet(workbook, "Ligand scores", 154)
        mappings = read_sheet(workbook, "Ligand classes", 14)
        classes = read_sheet(workbook, "Class scores", 33)
    finally:
        workbook.close()

    z_scores, assignments = check_ligand_tables(runs, scores, mappings)
    class_rows = check_class_tables(classes, z_scores, assignments)
    r = check_supplementary_tables(class_rows, z_scores)
    print("Raw runs: 1,232 recorded scores (880 original, 352 added)")
    print("Ligand scores: 154 minima and 154 population z-scores match")
    print("Class scores: 132 means, 33 recorded calls and 33 margins match")
    print("Selected rows: 15 of 15 agree with eligibility, margin ordering and benchmark calls")
    print("Table S7: 11 matched-ligand comparisons agree at published precision")
    print(f"Figure S3b: full-precision Pearson r = {r:.6f} (published -0.54)")
    print("Table S8: 4 candidates x 2 panels agree on the call, rounded margin and selection")


if __name__ == "__main__":
    main()
