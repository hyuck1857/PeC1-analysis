# Recorded docking scores and repeatability

This directory reproduces the study archive's calculations on completed Vina and Amber output. The two published examples contain original score lines. The seed table contains all 1,232 recorded scores from the original panel and its extension.

## Read the score examples

From the repository root, run with Python's standard library.

```sh
python docking/read_scores.py vina docking/data/vina_score_records.txt
python docking/read_scores.py mmgbsa docking/data/mmgbsa_difference.txt
```

The Vina example contains all nine `REMARK VINA RESULT` lines for the first archived original-panel row. The routine takes the first recorded mode, which gives **−5.047 kcal/mol**. The Amber example contains the `Differences (Complex - Receptor - Ligand)` section for the first successful row of the archived trimer MM/GBSA results. Its total is **−37.2450 kcal/mol**. The van der Waals, electrostatic, generalized Born and surface terms are −11.7383, −101.4928, 80.1339 and −4.1478 kcal/mol, respectively.

The original readers use the decimal notation found in the archived files. A partially populated Amber differences section yields `null` for missing terms. The command-line interface reports an error when it finds no readable score.

## Check variation across eight seeds

```sh
python docking/check_seed_variation.py
```

The calculation uses pandas, tested with version 2.3.3. Add `--output reproduced.csv` to save the complete recalculated table. `data/recorded_seed_scores.csv` has one row per original receptor–ligand pair, identified by a neutral `series` ID. Columns `seed_1` through `seed_8` preserve each recorded value and its seed. IDs were assigned within each collection after sorting by the eight scores in seed order. The original 110 series are shared across both panels. The extension adds 44 series.

| Archived panel | Series | Scores | Median within-series range, kcal/mol |
| --- | ---: | ---: | ---: |
| `original10` | 110 | 880 | 0.096 |
| `control14` | 154 | 1,232 | 0.092 |

For each series the original function computes count, best (minimum), median, mean, sample standard deviation, minimum, maximum, quartiles, interquartile range, range, and best minus median. The command compares every value with `data/archived_seed_statistics.csv` at absolute tolerance `1e-12`. The reference table carries the original archive values with the same neutral pair IDs. These statistics measure variability across computational seeds.

## Original code and source records

Source paths below are relative to the archived `SNU_phage_host_receptor_modeling` directory. The Python excerpts were checked against the original text with newline normalization. Imports, command-line interfaces and verification code were added for this deposit.

| Released routine | Original file and inclusive lines | Text preserved |
| --- | --- | --- |
| `read_scores.py` · `read_vina_result` | `15_blind_reanalysis/ultrablind_homologyfree/s14b_ultrablind_esmfold_vina.py`, 104–108 | The complete result-reading tail of `vina`, under a new function name. |
| `read_scores.py` · `parse_delta` | `scripts/s19_trimer_mmgbsa.py`, 157–164 | The complete function. The same definition occurs at lines 154–161 of `scripts/s20_trimer_mmgbsa_truncated.py`. |
| `check_seed_variation.py` · `pair_seed_statistics` | `CODEX/scripts/robustness_analysis.py`, 74–93 | The complete function, including pandas grouping and quantile calculations. The new caller uses the neutral series ID as one grouping key and a constant as the second. |

The original and additional seed records are from `out/raw_all_seeds.csv` and `out_with_control/raw_all_seeds.csv`, respectively, under `15_blind_reanalysis/ultrablind_homologyfree/`. Historical summary values come from `CODEX/outputs/seed_pair_statistics.csv`. For release verification, the Vina reader was also applied to the archived pose output for every seed row. All 880 original and 352 additional first-mode scores matched the corresponding CSV records.

The Amber example and comparison table come from `08_docking_and_binding_energy/definitive_af3_trimer/mmgbsa_trimer/`. The five parsed energy terms matched the archived `trimer_mmgbsa_results.csv` for all 18 rows marked `OK`. Six other rows are marked `FAILED` and have no reported total. Verification used existing reports throughout.

### Original file fingerprints

| Source | SHA-256 of original file |
| --- | --- |
| `s14b_ultrablind_esmfold_vina.py` | `803016e16f5be55fcad9cf72e194d69a613bd1fb552e81d992af35402c90397c` |
| `s19_trimer_mmgbsa.py` | `3a05f63acc31ffd0455afc146f3e14d115a8ee82e05610f384a1221b648fdbb9` |
| `s20_trimer_mmgbsa_truncated.py` | `ebf1fc87987d189e2652723f7bf5fa6c3d567f7016d4fb2e48216f6c545e1cc5` |
| `robustness_analysis.py` | `71435cf5c22d728b46a1c79bbc823f2f9c8386cc18bd34d5948764b396376388` |
| `out/raw_all_seeds.csv` | `0b201638b3aa550b02b81a8539130dd3b7845addbb51d2ecd212c80889f11aaf` |
| `out_with_control/raw_all_seeds.csv` | `e8a096d775725b32067bd7684570411a09ff31f854e29321402c457ec5cfc366` |
| `seed_pair_statistics.csv` | `efb803777d664b8c23a7e28ed4f703cd85295809065172ed696b09fed5894e6e` |
| Full file used for the Vina example | `ae00ee0b31577998c63e403366eb377a9f63adf76e3e63ad8b789fe96ead7fe1` |
| Full file used for the Amber example | `14326c7d39376924b422b50bd7cb9d84fb335c22ec4beba5380c3fe9f4a75ef8` |
| `trimer_mmgbsa_results.csv` | `34f818c90b453144de2cd1de8b5c56a404e9ea28e3ad26467989c779bc03470e` |
