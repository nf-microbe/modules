#!/usr/bin/env python

import argparse
import gzip
import os
import sys

import pandas as pd
from Bio import SeqIO


def parse_args(args=None):
    description = "Filter sequences based on classification, composition, and completeness thresholds."
    epilog = """
    Example usage:
    python filtersequences.py \
        --input_fasta sequences.fasta.gz \
        --input_proteins proteins.faa.gz \
        --completeness_data completeness_data.tsv.gz \
        --seqs_to_keep seqs_to_keep.tsv \
        --output_fasta filtered.fasta \
        --output_proteins filtered_proteins.faa
    """

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-if",
        "--input_fasta",
        help="Path to input FASTA file (gzipped) containing unfiltered sequences.",
    )
    parser.add_argument(
        "-ip",
        "--input_proteins",
        help="Path to input FAA file (gzipped) containing unfiltered protein sequences.",
    )
    parser.add_argument(
        "-c",
        "--completeness_data",
        help="Path to TSV file (gzipped) containing completeness data.",
    )
    parser.add_argument(
        "-s",
        "--seqs_to_keep",
        help="Path to TSV file containing sequence IDs to keep.",
    )
    parser.add_argument(
        "-of",
        "--output_fasta",
        help="Path to output FASTA file containing filtered sequences.",
    )
    parser.add_argument(
        "-op",
        "--output_proteins",
        help="Path to output FAA file containing filtered protein sequences.",
    )
    return parser.parse_args(args)


def filter_sequences(
    input_fasta,
    input_proteins,
    completeness_data,
    seqs_to_keep,
    output_fasta,
    output_proteins,
):
    """
    Filter sequences based on classification, composition, and completeness thresholds.

    Args:
        input_fasta         : Path to input FASTA file (gzipped) containing unfiltered sequences.
        input_proteins      : Path to input FAA file (gzipped) containing unfiltered protein sequences.
        completeness_data   : Path to TSV file (gzipped) containing completeness data.
        seqs_to_keep        : Path to TSV file containing sequence IDs to keep.
        output_fasta        : Path to output FASTA file containing filtered sequences.
        output_proteins     : Path to output FAA file containing filtered protein sequences.

    Returns:
        Outputs FASTA file and FAA file with filtered sequences/proteins
    """
    if os.path.exists(seqs_to_keep):
        # Read in sequence IDs to keep
        override_seqs_df = pd.read_csv(seqs_to_keep, sep="\t", names=["seqid"], header=None)
        override_seqs_set = set(override_seqs_df["seqid"].tolist())
    else:
        override_seqs_set = set()

    # Read in completeness data
    completeness_data_df = pd.read_csv(
        completeness_data,
        sep="\t",
        header=0,
        usecols=["contig_id", "completeness_status", "composition_status", "virus_classification"],
    )

    # identify sequences passing all filters
    passing_set = set(
        completeness_data_df[
            (completeness_data_df["completeness_status"] == "pass")
            & (completeness_data_df["composition_status"] == "pass")
            & (completeness_data_df["virus_classification"] == "viral")
        ]["contig_id"].to_list()
    )
    seqs_to_keep_set = override_seqs_set.union(passing_set)

    # Write out filtered sequences
    sequences_to_write = []
    fasta = gzip.open(input_fasta, "rt") if input_fasta.split(".")[-1] == "gz" else open(input_fasta)
    for record in SeqIO.parse(fasta, "fasta"):
        if record.id in seqs_to_keep_set:
            sequences_to_write.append(record)
    SeqIO.write(sequences_to_write, output_fasta, "fasta")

    # Write out filtered protein sequences
    proteins_to_write = []
    if os.path.exists(input_proteins):
        proteins = gzip.open(input_proteins, "rt") if input_fasta.split(".")[-1] == "gz" else open(input_proteins)
        for record in SeqIO.parse(proteins, "fasta"):
            if record.id.rpartition("_")[0] in seqs_to_keep_set:
                proteins_to_write.append(record)
        SeqIO.write(proteins_to_write, output_proteins, "fasta")
    else:
        open(output_proteins, "a").close()


def main(args=None):
    args = parse_args(args)
    filter_sequences(
        args.input_fasta,
        args.input_proteins,
        args.completeness_data,
        args.seqs_to_keep,
        args.output_fasta,
        args.output_proteins,
    )


if __name__ == "__main__":
    sys.exit(main())
