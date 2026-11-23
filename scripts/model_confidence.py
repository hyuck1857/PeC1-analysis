"""Recalculate the correlation in Figure S3b from the values in Table S7.

The averaging and Pearson correlation follow fig3 in the study figure script
make_figs_a.py. The CSV reader and command-line output were added for review.
"""

import csv
from pathlib import Path

import numpy as np


DATA = Path(__file__).resolve().parents[1] / "data" / "model_confidence.csv"


def main():
    with DATA.open(encoding="utf-8", newline="") as source:
        rows = list(csv.DictReader(source))
    plddt = np.array([float(row["local_plddt"]) for row in rows])
    meanz = np.array([(float(row["z_cellotriose"]) + float(row["z_galactan"])) / 2
                      for row in rows])
    r = np.corrcoef(plddt, meanz)[0, 1]
    print(f"Candidates: {len(rows)}")
    print(f"Pearson r (local pLDDT vs mean matched-ligand z-score): {r:.2f}")
    if len(rows) != 11 or round(float(r), 2) != -0.54:
        raise ValueError("The correlation differs from the value reported in Figure S3b.")


if __name__ == "__main__":
    main()
