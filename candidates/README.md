# Benchmark nominations (Table S1)

Run `python candidates/check_table_s1.py` from the repository root.

`archived_cohort.csv` lists all 859 benchmark CDS records. For each working phage name and locus tag, the final column records membership in the archived 85-row output. The identifiers come from `02_SNU_phage_ORF_annotation.csv` (SHA-256 `d3527f6fd31e2c1630254408184f0a3fdb131d8104fe72d735400ec807de8750`). The checker verifies the 85 recorded memberships against the candidate table.

`archived_candidate_priorities.csv` contains five columns from the original 85-row output, `03_SNU_phage_RBP_candidates.csv`: working phage name, locus tag, recorded candidate priority, protein length and the working product annotation. These records precede the selection of a representative for each phage. The priority values in this file are annotation categories for candidate proteins. The information tiers in manuscript Table 2 refer to a separate comparison.

The script uses the sorting and counting operations from `s01_inventory_qc_annotation.py`. Within each phage, that code orders candidates by ascending priority number and then by descending protein length. The new comparison reads the first candidate in each group and verifies its locus tag and length against the five rows of final Table S1. The working label phiEaP-33 is mapped to manuscript label phiEaP-21 for that comparison. Priority assignment from genome annotations takes place before the input to this script.

Source fingerprints (SHA-256) are `50010ff0b0cc50872a4aa033deaf9e14fe2f4941010335e8bc55b82886202765` for the historical candidate table and `df3b0e626f9bae25b7ada92cf1f98c8e5e71cd4b61cb44cce66a30e071bbd405` for the original source script.
