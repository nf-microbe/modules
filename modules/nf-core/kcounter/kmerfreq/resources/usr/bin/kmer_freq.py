#!/usr/bin/env python

import argparse
import gzip
import multiprocessing as mp
import sys

import kcounter
import numpy as np
from Bio import SeqIO


def parse_args(args=None):
    description = "Calculate kmer frequency for nucleotide sequences."
    epilog = "Example usage: python kmer_freq.py -i sequences.fasta.gz -o kcounter.tsv"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to input FASTA to calculate kmer frequency from",
    )
    parser.add_argument(
        "-t",
        "--threads",
        help="Number of threads to use in multiprocessing",
    )
    parser.add_argument(
        "-k",
        "--kmer_len",
        help="K-mer length to use when calculating frequencies",
        default=21,
        type=int,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output TSV file kmer frequencies.",
    )
    return parser.parse_args(args)


def get_average_kmer_freq(genome, kmer_len):
    counts = list(kcounter.count_kmers(str(genome.seq), kmer_len, canonical_kmers=True).values())
    return round(np.mean(counts), 2)


def main(args=None):
    args = parse_args(args)

    # read in genomes
    genomes = []
    input_gunzipped = gzip.open(args.input, "rt") if args.input.split(".")[-1] == "gz" else open(args.input)
    for record in SeqIO.parse(input_gunzipped, "fasta"):
        genomes.append(record)

    with mp.Pool(int(args.threads)) as pool:
        kmer_freq_list = pool.starmap(get_average_kmer_freq, zip(genomes, [int(args.kmer_len)] * len(genomes)))

    with open(args.output, "w") as fout:
        fout.write("contig_id\tkmer_freq\n")
        for genome, kmer_freq in zip(genomes, kmer_freq_list):
            fout.write(f"{genome.id}\t{kmer_freq}\n")


if __name__ == "__main__":
    sys.exit(main())
