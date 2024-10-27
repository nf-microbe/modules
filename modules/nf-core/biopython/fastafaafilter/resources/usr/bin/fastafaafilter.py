#!/usr/bin/env python

import argparse
import gzip
import sys

from Bio import SeqIO


def parse_args(args=None):
    description = "Filter FastA and FAA sequences based on length."
    epilog = """
    Example usage:
    python fastafaafilter.py \
        --input_fasta sequences.fasta.gz \
        --input_faa proteins.faa.gz \
        --prefix test \
        --min_fasta_len 1000
    """

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "--input_fasta",
        help="Path to input FASTA file (gzipped) containing unfiltered sequences.",
    )
    parser.add_argument(
        "--input_faa",
        help="Path to input FAA file (gzipped) containing unfiltered protein sequences.",
    )
    parser.add_argument(
        "--prefix",
        help="Prefix to use when naming output files.",
    )
    parser.add_argument(
        "--fasta_min_len",
        help="Minimum length to use when filtering FastA sequences.",
    )
    return parser.parse_args(args)


def filter_sequences(input_fasta, input_faa, prefix, fasta_min_len):
    """
    Filter sequences based on classification, composition, and completeness thresholds.

    Args:
        input_fasta         : Path to input FASTA file (gzipped) containing unfiltered sequences.
        input_faa           : Path to input FAA file (gzipped) containing unfiltered protein sequences.
        prefix              : File prefix to use when naming output files.
        fasta_min_len       : Minimum length of FastA sequence to retain.

    Returns:
        Outputs FASTA file and FAA file with filtered sequences/proteins
    """
    min_len = int(fasta_min_len)
    # Write out filtered sequences
    sequences_to_write = []
    passing_record_ids = set()
    fasta = gzip.open(input_fasta, "rt") if input_fasta.split(".")[-1] == "gz" else open(input_fasta)
    for record in SeqIO.parse(fasta, "fasta"):
        if len(record.seq) > min_len:
            sequences_to_write.append(record)
            passing_record_ids.add(record.id)
    SeqIO.write(sequences_to_write, prefix + ".fasta", "fasta")

    # Write out filtered protein sequences
    proteins_to_write = []
    proteins = gzip.open(input_faa, "rt") if input_fasta.split(".")[-1] == "gz" else open(input_faa)
    for record in SeqIO.parse(proteins, "fasta"):
        if record.id.rpartition("_")[0] in passing_record_ids:
            proteins_to_write.append(record)
    SeqIO.write(proteins_to_write, prefix + ".faa", "fasta")


def main(args=None):
    args = parse_args(args)
    filter_sequences(args.input_fasta, args.input_faa, args.prefix, args.fasta_min_len)


if __name__ == "__main__":
    sys.exit(main())
