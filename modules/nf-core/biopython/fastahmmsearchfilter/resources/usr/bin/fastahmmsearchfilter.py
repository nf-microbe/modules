#!/usr/bin/env python

import argparse
import gzip
import sys

from Bio import SeqIO


def parse_args(args=None):
    description = "Filter FastA sequences based on HMM hits."
    epilog = """
    Example usage:
    python fastahmmsearchfilter.py \
        --input_fasta sequences.fasta.gz \
        --hmm_tbl hmmsearch.tbl \
        --prefix test
    """

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "--input_fasta",
        help="Path to input FASTA file (gzipped) containing unfiltered sequences.",
    )
    parser.add_argument(
        "--hmm_tbl",
        help="Path to input HMMsearch's tblout file (gzipped) containing HMM hits.",
    )
    parser.add_argument(
        "--prefix",
        help="Prefix to use when naming output files.",
    )

    return parser.parse_args(args)


def filter_sequences(input_fasta, hmm_tbl, prefix):
    """
    Filter sequences based on classification, composition, and completeness thresholds.

    Args:
        input_fasta         : Path to input FASTA file (gzipped) containing unfiltered sequences.
        hmm_tbl             : Path to input HMMsearch's tblout file (gzipped) containing HMM hits.
        prefix              : File prefix to use when naming output files.

    Returns:
        Outputs FASTA file with only sequences containing HMM hits
    """
    # Identify sequences with HMM hits
    seqs_w_hmm_set = set()
    tbl_path_open = gzip.open(hmm_tbl, "rt") if hmm_tbl.split(".")[-1] == "gz" else open(hmm_tbl)
    for line in tbl_path_open:
        if line[0] == "#":
            continue
        parts = line.split()
        seqs_w_hmm_set.add(parts[0].rpartition('_')[0])

    # Write out filtered sequences
    sequences_to_write = []
    fasta = gzip.open(input_fasta, "rt") if input_fasta.split(".")[-1] == "gz" else open(input_fasta)
    for record in SeqIO.parse(fasta, "fasta"):
        if record.id in seqs_w_hmm_set:
            sequences_to_write.append(record)
    SeqIO.write(sequences_to_write, prefix + ".fasta", "fasta")


def main(args=None):
    args = parse_args(args)
    filter_sequences(args.input_fasta, args.hmm_tbl, args.prefix)


if __name__ == "__main__":
    sys.exit(main())
