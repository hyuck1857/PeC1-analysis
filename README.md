# PeC1 analysis

Companion code and numerical data for *Computational prediction and experimental validation of host requirements for receptor recognition by Pectobacterium phage PeC1*.

This repository covers genome quality control, assembly depth, benchmark scoring and two model-confidence analyses. Reviewers can read the files on GitHub or download them using **Code → Download ZIP**, without signing in.

## Reproduce the reported summaries

Run these commands from the repository directory. The benchmark and interface checks use only the Python standard library. The correlation calculation uses NumPy. The genome QC utility uses Biopython.

```sh
python -m pip install -r requirements.txt
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

The deposit starts from recorded calls and model scores for the benchmark checks. Upstream candidate nomination, structure generation, docking runs and the full receptor-prediction workflow are outside its scope.
