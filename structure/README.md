# Structure prediction inputs, model files and confidence

[original_helpers.py](original_helpers.py) contains seven functions copied verbatim from the study archive. [model_files.py](model_files.py) provides a command-line interface added for this deposit. Run the examples from the repository root with the dependencies in `requirements.txt`.

## Reproduce the Figure S3a calibration

```sh
python structure/model_files.py replay
```

Expected output is `Matched all 31 job scores from 155 model summaries.` Each job contributed five AlphaFold Server summaries. The original `best_model_metrics` function selects the highest `ranking_score`, retaining the first filename in sorted order when scores tie. The replay checks the resulting ipTM values against `data/interface_confidence.csv` within each group, including repeated observations.

[calibration_confidences.zip](calibration_confidences.zip) contains the original JSON file contents from the 31 jobs listed in the archived `07_calibration_master.csv`. The deposit relabels their paths with the four public group names and a job number. Job numbers follow the master table's order within each group. Each file contains numeric confidence fields from one model. Its final filename digit is the model index.

## Work with a single model

| Command | Input and output |
| --- | --- |
| `python structure/model_files.py af-json protein.fasta --name example --copies 1` | Reads a single protein FASTA record and prints a one-job JSON array in the archived AlphaFold Server format. The caller supplies the copy count. |
| `python structure/model_files.py convert model.cif model.pdb` | Converts the chosen mmCIF to PDB and reports the converter used. The historical function accepts output files larger than 1,000 bytes. The command writes the PDB after conversion succeeds. |
| `python structure/model_files.py metrics completed_job` | Reads the summary JSONs in one completed job directory and prints the best-ranked model's ipTM, pTM, ranking score, interface ipTM and number of models. |
| `python structure/model_files.py plddt model.pdb` | Reads C-alpha B-factor fields from ESMFold `ATOM` records and reports the arithmetic mean to one decimal. The original code multiplies the mean by 100 when it is at most 1.5 and returns 0.0 for an empty C-alpha selection. |

The original `esmfold(sequence, outpdb)` function sends one protein to the historical ESM Atlas endpoint when invoked with an uncached output path. The study code reuses an existing file larger than 1,000 bytes. The public command-line examples run offline. An offline test checks the request body, response write and cache behavior with a mocked reply. Live service availability was untested.

```sh
python structure/test_original_helpers.py
```

## Archived sources

Line numbers are one-based and inclusive. Only the module-level imports and opening docstring in the deposited helper module were supplied for this release.

| Archived source and original line numbers | Function and use in the study |
| --- | --- |
| `SNU_phage_host_receptor_modeling/scripts/s13_alphafold_server_json.py`, 36–37 and 43–45 | `safe_name`, `job`. Format a job name and build the JSON object for one protein. Enclosing the returned object in a list reproduces the saved single-job format. |
| `SNU_phage_host_receptor_modeling/15_blind_reanalysis/ultrablind_homologyfree/s14b_ultrablind_esmfold_vina.py`, 58–69 and 71–75 | `esmfold`, `mean_plddt`. The original Tier 3 single-protein request and the subsequent confidence calculation across the submitted model's C-alpha records. |
| `SNU_phage_host_receptor_modeling/scripts/s17_import_af_trimers.py`, 50–53 and 55–80 | `rank_of`, `cif_to_pdb`. The importer chose the lowest server model index from the CIF filenames. It tried Gemmi, Biopython and Open Babel in order. The study chronicle records that Biopython performed the eight primary conversions. |
| `C1_interaction_boost_analysis/scripts/06_parse_calibration_results.py`, 51–64 | `best_model_metrics`. Select one model per completed calibration job using the numerical ranking score in its summary JSON. |

## Checks against the study archive

- The job formatter reproduced all 57 saved single-job JSON objects from their filenames, stored protein inputs and copy counts.
- Each of the eight primary mmCIFs matched a `model_0.cif` in the archived server downloads byte for byte. With Biopython 1.81, conversion reproduced all eight archived PDBs after normalizing line endings. All 83,973 atom B-factors matched the source mmCIFs at PDB precision.
- For all 11 archived ESMFold PDBs, `mean_plddt` reproduced the archived metadata and the local pLDDT values reported in Table S7. The four records archived under phiEaP-33 appear under phiEaP-21 in the manuscript.
- For each of the 31 calibration jobs, all five fields returned by `best_model_metrics` matched the archived calibration master table. All 155 JSON files in the deposit matched their corresponding archived file contents byte for byte.
