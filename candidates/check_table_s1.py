"""Replay the archived candidate ordering and check the nominations in Table S1."""

from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent


def main():
    records = pd.read_csv(HERE / "archived_candidate_priorities.csv")
    if len(records) != 85 or records.duplicated(["phage", "locus_tag"]).any():
        raise ValueError("Expected the 85 distinct records from the archived candidate list.")
    cohort = pd.read_csv(HERE / "archived_cohort.csv")
    keys = ["phage", "locus_tag"]
    if len(cohort) != 859 or cohort.duplicated(keys).any():
        raise ValueError("The original benchmark contains 859 distinct CDS records.")
    if not cohort["in_archived_candidate_list"].isin([0, 1]).all():
        raise ValueError("Invalid archived candidate-list membership flag.")
    included = cohort.loc[cohort["in_archived_candidate_list"] == 1, keys]
    pd.testing.assert_frame_equal(included.sort_values(keys).reset_index(drop=True),
                                  records[keys].sort_values(keys).reset_index(drop=True))
    rbp_rows = records.to_dict(orient="records")
    # Sorting and counting operations from s01_inventory_qc_annotation.py.
    rbp = pd.DataFrame(rbp_rows).sort_values(["phage","candidate_tier","length_aa"],
                                             ascending=[True, True, False]).reset_index(drop=True)
    summary = rbp.groupby(["phage","candidate_tier"]).size().unstack(fill_value=0)
    print("Archived counts by candidate priority")
    print(summary.to_string())

    # Compare the first entry in each phage group with final Table S1.
    selected = rbp.drop_duplicates("phage")[["phage", "locus_tag", "length_aa"]].copy()
    selected["phage"] = selected["phage"].replace({"phiEaP-33": "phiEaP-21"})
    selected["locus_tag"] = selected["locus_tag"].str.replace(
        "phiEaP-33_", "phiEaP-21_", regex=False
    )
    expected = pd.read_csv(HERE / "table_s1_nominations.csv")
    pd.testing.assert_frame_equal(
        selected.sort_values("phage").reset_index(drop=True),
        expected.sort_values("phage").reset_index(drop=True),
    )
    print("\nAll five locus tags and lengths agree with Table S1.")


if __name__ == "__main__":
    main()
