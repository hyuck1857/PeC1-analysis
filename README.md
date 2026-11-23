# PeC1 analysis

Companion code and numerical data for *Computational prediction and experimental validation of host requirements for receptor recognition by Pectobacterium phage PeC1*.

This repository contains study code and archived data for the benchmark nominations, structure models, docking scores and PeC1 assembly statistics. Reviewers can [download the review-v2 snapshot](https://github.com/hyuck1857/PeC1-analysis/archive/refs/tags/review-v2.zip) without signing in.

## Trace a result back to its inputs

| Manuscript result | Files and scope |
| --- | --- |
| Five nominations in Table S1 | [Candidate records and ordering code](candidates/). The cohort lists all 859 CDS records and identifies the 85 originally nominated candidates. Applying the historical sort to their recorded priorities and lengths recovers all five Table S1 entries. The initial priority assignment is upstream of this check. |
| Structure prediction and confidence | [Original structure helpers and model summaries](structure/). The functions include single-protein AlphaFold Server job specification, the archived ESMFold request, structure-file conversion and model-confidence extraction. Numerical summaries for all five models in each of 31 calibration jobs allow the reviewer to repeat top-model selection offline. |
| Docking repeatability | [Docking output readers and seed statistics](docking/). The archived code excerpts read Vina and Amber outputs and summarize the variation across random seeds. |
| Tier 3 scores and calls | [Final Supplementary Data S1](data/Supplementary_data_S1_v16_20261102.xlsx) and [workbook checker](scripts/check_supplementary_s1.py). The checker starts with all 1,232 recorded docking scores, verifies 154 minima and standardized scores, and reconstructs 33 candidate-level calls and all 15 phage-level selections across the three panel conditions. |

## Reproduce the reported summaries

Run these commands from the repository directory. The review checks were tested with Python 3.9.0 and the versions in `requirements.txt`.

```sh
python -m pip install -r requirements.txt
python candidates/check_table_s1.py
python structure/model_files.py replay
python scripts/check_supplementary_s1.py
python docking/check_seed_variation.py
python scripts/check_benchmark.py
python scripts/check_interface_confidence.py
python scripts/model_confidence.py
python scripts/summarize_coverage.py --frequencies data/assembly_depth_frequency.csv
```

| Manuscript item | Input and calculation | Expected result |
| --- | --- | --- |
| Table 2 | Score the five sets of recorded calls in `data/benchmark_calls.csv` against the accepted classes for the three graded phages. Sensitivity calls come from the selected rows in Supplementary Data S1. | Tier 2 recovers 3/3. Tier 3 recovers 2/3 with the original panel, 2/3 in the expanded control and 1/3 in the expanded sensitivity analysis. The best fixed class recovers 2/3. |
| Figure S3a and Table S5 | Summarize the highest-ranked model for each of the 31 calibration jobs in `data/interface_confidence.csv`. Scores are grouped and retain repeated observations. | Positive controls span 0.13–0.23 and negative controls span 0.10–0.14. One positive control falls within the negative-control range. |
| Figure S3b and Table S7 | Correlate local pLDDT with the mean of the two matched-ligand z-scores in `data/model_confidence.csv`, transcribed at the precision reported in Table S7. | Pearson r = −0.54 across 11 candidates. |
| PeC1 assembly depth (Table 3) | Compute weighted summary statistics from the number of positions at each observed depth in `data/assembly_depth_frequency.csv`. This histogram was made from all 40,178 positions of the assembly contig used for read mapping. | Mean 17,951.7×, median 17,630×, range 4,375–55,303× and zero uncovered positions. The final genome after collapse of the terminal overlap is 40,037 bp. |

## Genome quality control

Supply a single-record nucleotide FASTA and its matching annotated GenBank record.

```sh
python scripts/genome_qc.py --fasta genome.fasta --genbank genome.gbk
```

The script writes one CSV row with sequence length, GC percentage, CDS count, summed CDS span as a percentage of genome length and ambiguous-base counts. The calculations were adapted from `s01_inventory_qc_annotation.py` and checked against the archived QC output for all five working inputs. One historical input labelled phiEaP-33 corresponds to the manuscript entry phiEaP-21. The archived record gives 84,662 bp and 43.82% GC. Final Table 1 lists 84,541 bp and 43.83% GC. Both records report 117 CDS and 88.07% coding density. The archived input and the final table therefore give different sequence-length and GC values.

## Code provenance and scope

The genome QC, ipTM group summaries and pLDDT correlation preserve calculations from the study scripts `s01_inventory_qc_annotation.py`, `06_parse_calibration_results.py` and `make_figs_a.py`, respectively. They have standalone entry points for reviewers. The Table 2 checker was written for this deposit to apply the final scoring rule to the reported calls. It accepts either colanic acid or LPS for phiPccP-2 and excludes the two ungraded phages from the denominator.

The component READMEs identify original functions and later verification code. For candidate selection, the reproducible starting point is the archived list with assigned priorities. For docking, it is the archived outputs and individual scores. The genome-to-candidate prioritization pipeline and the workflow that prepares targets and executes new docking runs are not deposited. The structure section identifies which operations were verified offline and which originally used external services.
