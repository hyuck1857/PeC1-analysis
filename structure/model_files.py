"""Run the deposited study functions on single files and the recorded calibration."""

import argparse
import csv
import json
import os
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import Path as ZipPath, ZipFile

from original_helpers import best_model_metrics, cif_to_pdb, job, mean_plddt


HERE = Path(__file__).resolve().parent


def replay_calibration():
    with (HERE.parent / "data/interface_confidence.csv").open(newline="") as source:
        expected = Counter((row["group"], float(row["iptm"])) for row in csv.DictReader(source))

    observed = Counter()
    model_count = 0
    with ZipFile(HERE / "calibration_confidences.zip") as archive:
        for group in ZipPath(archive).iterdir():
            for run in group.iterdir():
                with TemporaryDirectory() as tmp:
                    for index in range(5):
                        name = f"model_summary_confidences_{index}.json"
                        (Path(tmp) / name).write_bytes((run / name).read_bytes())
                    metrics = best_model_metrics(tmp)
                observed[group.name, metrics["iptm"]] += 1
                model_count += metrics["n_models"]
        if len(archive.namelist()) != model_count:
            raise ValueError("Archive contains unexpected files.")
    if observed != expected:
        raise ValueError("Recalculated group scores differ from data/interface_confidence.csv.")
    print(f"Matched all {sum(observed.values())} job scores from {model_count} model summaries.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("replay", help="Recalculate the 31 deposited calibration scores")
    af = commands.add_parser("af-json", help="Print the archived AlphaFold Server JSON format for one protein")
    af.add_argument("fasta", type=Path)
    af.add_argument("--name", required=True)
    af.add_argument("--copies", type=int, default=1)
    convert = commands.add_parser("convert", help="Convert one model mmCIF to a PDB")
    convert.add_argument("cif", type=Path)
    convert.add_argument("pdb", type=Path)
    confidence = commands.add_parser("metrics", help="Read confidence summaries for one completed AF3 job")
    confidence.add_argument("job_directory", type=Path)
    plddt = commands.add_parser("plddt", help="Calculate mean C-alpha pLDDT for an ESMFold PDB")
    plddt.add_argument("pdb", type=Path)
    args = parser.parse_args()

    if args.command == "replay":
        replay_calibration()
    elif args.command == "af-json":
        from Bio import SeqIO
        if args.copies < 1:
            parser.error("--copies must be positive.")
        sequence = str(SeqIO.read(args.fasta, "fasta").seq)
        print(json.dumps([job(args.name, sequence, args.copies)], indent=2))
    elif args.command == "convert":
        if args.cif.resolve() == args.pdb.resolve():
            parser.error("Choose a different path for the PDB output.")
        with TemporaryDirectory(dir=args.pdb.parent) as tmp:
            result = Path(tmp) / "converted.pdb"
            method = cif_to_pdb(str(args.cif), str(result))
            if method is None:
                parser.error("The archived converter could not produce a PDB larger than 1,000 bytes.")
            os.replace(result, args.pdb)
        print(f"Converted with {method}.")
    elif args.command == "metrics":
        metrics = best_model_metrics(str(args.job_directory))
        if metrics is None:
            parser.error("The job directory contains no matching summary confidence files.")
        print(json.dumps(metrics, indent=2))
    elif args.command == "plddt":
        print(f"{mean_plddt(str(args.pdb)):.1f}")


if __name__ == "__main__":
    main()
