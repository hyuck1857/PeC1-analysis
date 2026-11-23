"""Summarize observed read depth for one reference from an existing table.

This is a verification script for the archived summary of the 40,178 bp
pre-collapse assembly contig. The 40,037 bp final genome has a different
coordinate system. The reported mapping depth applies to the former.
"""

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path


def read_position_depths(handle):
    """Read a complete, headerless reference, position, depth table."""
    counts = Counter()
    reference = None
    for expected_position, line in enumerate(handle, start=1):
        fields = line.split()
        if len(fields) != 3:
            raise ValueError(f"Expected reference, position and depth on line {expected_position}.")
        name, position, depth = fields
        if reference is None:
            reference = name
        if name != reference or int(position) != expected_position:
            raise ValueError("Expected one reference with consecutive positions starting at 1.")
        depth = int(depth)
        if depth < 0:
            raise ValueError("Depth cannot be negative.")
        counts[depth] += 1
    return counts


def read_frequencies(handle):
    """Read a CSV with the number of reference positions at each depth."""
    counts = Counter()
    reader = csv.DictReader(handle)
    if reader.fieldnames != ["depth", "positions"]:
        raise ValueError("The expected CSV header is depth,positions.")
    for row in reader:
        if None in row or None in row.values():
            raise ValueError(f"Expected two columns on line {reader.line_num}.")
        depth, positions = int(row["depth"]), int(row["positions"])
        if depth < 0 or positions < 1 or depth in counts:
            raise ValueError(f"Invalid or duplicate depth-frequency entry on line {reader.line_num}.")
        counts[depth] = positions
    return counts


def summarize(counts):
    positions = sum(counts.values())
    if positions == 0:
        raise ValueError("Input contains no observed positions.")

    # The two median ranks coincide when the number of positions is odd.
    left_rank, right_rank = (positions + 1) // 2, (positions + 2) // 2
    cumulative = 0
    for depth, frequency in sorted(counts.items()):
        if cumulative < left_rank <= cumulative + frequency:
            left_depth = depth
        cumulative += frequency
        if cumulative >= right_rank:
            median = (left_depth + depth) / 2
            break

    total_depth = sum(depth * frequency for depth, frequency in counts.items())
    return {
        "reference_positions": positions,
        "mean_depth": round(total_depth / positions, 1),
        "median_depth": median,
        "min_depth": min(counts),
        "max_depth": max(counts),
        "zero_depth_positions": counts[0],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--frequencies", type=Path, help="CSV with header depth,positions")
    group.add_argument("--depth-table", type=Path, help="Reference, position and depth, including zeros")
    args = parser.parse_args()
    try:
        with (args.frequencies or args.depth_table).open(encoding="utf-8-sig", newline="") as handle:
            counts = read_frequencies(handle) if args.frequencies else read_position_depths(handle)
        summary = summarize(counts)
    except (ValueError, OSError, csv.Error) as exc:
        parser.exit(1, f"Coverage summary failed: {exc}\n")
    writer = csv.DictWriter(sys.stdout, fieldnames=list(summary), lineterminator="\n")
    writer.writeheader()
    writer.writerow(summary)


if __name__ == "__main__":
    main()
