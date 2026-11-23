"""Reproduce the archived within-pair score statistics for both seed panels."""

import argparse
from pathlib import Path

import pandas as pd


def pair_seed_statistics(raw: pd.DataFrame, panel: str) -> pd.DataFrame:
    grouped = raw.groupby(["receptor", "ligand"])["affinity"]
    stats = grouped.agg(
        n_seeds="count",
        best="min",
        median="median",
        mean="mean",
        std="std",
        minimum="min",
        maximum="max",
    ).reset_index()
    q = grouped.quantile([0.25, 0.75]).unstack()
    q.columns = ["q25", "q75"]
    q = q.reset_index()
    stats = stats.merge(q, on=["receptor", "ligand"], how="left")
    stats["iqr"] = stats["q75"] - stats["q25"]
    stats["range"] = stats["maximum"] - stats["minimum"]
    stats["best_minus_median"] = stats["best"] - stats["median"]
    stats.insert(0, "panel", panel)
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write the recomputed table as CSV")
    args = parser.parse_args()
    data = Path(__file__).resolve().parent / "data"
    wide = pd.read_csv(data / "recorded_seed_scores.csv")
    seed_columns = [f"seed_{seed}" for seed in range(1, 9)]
    panels = {"original10": wide[wide["collection"] == "original"], "control14": wide}
    records = []
    for panel, frame in panels.items():
        raw = frame.melt(id_vars="series", value_vars=seed_columns, value_name="affinity")
        # Each neutral series ID represents one original receptor-ligand pair.
        raw = raw.rename(columns={"series": "receptor"}).assign(ligand="recorded")
        stats = pair_seed_statistics(raw, panel)
        records.append(stats.rename(columns={"receptor": "series"}).drop(columns="ligand"))
    calculated = pd.concat(records, ignore_index=True)
    expected = pd.read_csv(data / "archived_seed_statistics.csv")
    pd.testing.assert_frame_equal(calculated, expected, rtol=0, atol=1e-12)
    if args.output:
        calculated.to_csv(args.output, index=False)
    for panel, group in calculated.groupby("panel", sort=False):
        print(f"{panel}: {len(group)} series, {group['n_seeds'].sum()} scores, "
              f"median within-series range {group['range'].median():.3f} kcal/mol")
    print(f"Verified all {len(calculated)} archived panel-series rows.")


if __name__ == "__main__":
    main()
