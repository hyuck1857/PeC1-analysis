"""Report genome QC for one matching pair of single-record FASTA and GenBank files.

Adapted from the benchmark workflow's s01_inventory_qc_annotation.py. The GC,
ambiguity and CDS-span calculations follow that script. The original workflow
read the coding-density denominator from its sample configuration. Here it is
read from the supplied FASTA, whose length matched the configuration for all
five archived benchmark FASTA files.
"""

import argparse
import csv
import sys

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


def genome_qc(fasta, genbank):
    record = SeqIO.read(fasta, "fasta")
    annotation = SeqIO.read(genbank, "genbank")
    sequence = str(record.seq).upper()
    genome_length = len(sequence)
    if not genome_length:
        raise ValueError("The FASTA record contains no sequence.")
    if sequence != str(annotation.seq).upper():
        raise ValueError("The supplied FASTA and GenBank sequences differ.")

    gc = round(100 * gc_fraction(record.seq), 2)
    ambiguous = sum(sequence.count(base) for base in "RYSWKMBDHVN")
    cds = [feature for feature in annotation.features if feature.type == "CDS"]
    total_cds_span = 0
    for feature in cds:
        if feature.location is None:
            raise ValueError("A CDS in the GenBank record has no location.")
        # GenBank coordinates are parsed as zero-based, half-open locations.
        # As in the original workflow, a joined feature contributes its full
        # outer span, and overlaps between separate CDS features count twice.
        total_cds_span += int(feature.location.end) - int(feature.location.start)

    return {
        "contig_id": record.id,
        "genome_length_bp": genome_length,
        "n_contigs": 1,
        "n_CDS": len(cds),
        "GC_percent": gc,
        "coding_density_percent": round(100 * total_cds_span / genome_length, 2),
        "N_bases": sequence.count("N"),
        "ambiguous_bases": ambiguous,
        "ambiguous_fraction": round(ambiguous / genome_length, 6),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--fasta", required=True, help="Single-record nucleotide FASTA file")
    parser.add_argument("--genbank", required=True, help="Matching annotated GenBank file")
    args = parser.parse_args()
    try:
        metrics = genome_qc(args.fasta, args.genbank)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"Genome QC failed: {exc}\n")
    writer = csv.DictWriter(sys.stdout, fieldnames=list(metrics), lineterminator="\n")
    writer.writeheader()
    writer.writerow(metrics)


if __name__ == "__main__":
    main()
