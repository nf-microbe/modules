#!/usr/bin/env python

import argparse
import gzip
import sys


def parse_args(args=None):
    description = "Multiply a sequence based on the average kmer abundance displayed in the header."
    epilog = "Example usage: python multiplier.py --input sequences.fasta --output test.fasta"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to input FASTA (unzipped) to be multiplied.",
    )
    parser.add_argument(
        "-a",
        "--abundance",
        help="Minimum contig abundance to include sequence in output FastA.",
        default=0,
        type=float,
    )
    parser.add_argument("-m", "--multiplier", help="Multiply sequences by kmer abundance", action="store_true")
    parser.add_argument(
        "-o",
        "--output",
        help="Output FastA files containing multiplied sequences.",
    )
    return parser.parse_args(args)


def parse_fasta(handle):
    header = next(handle).strip()[1:]
    id = header.split()[0]
    mult = float(header.split()[1].split(":")[2])
    seq = ""
    for line in handle:
        if line[0] == ">":
            yield header, id, mult, seq
            header = line.strip()[1:]
            id = header.split()[0]
            mult = float(header.split()[1].split(":")[2])
            seq = ""
        else:
            seq += line.rstrip().upper()
    yield header, id, mult, seq


def multiply_sequences(input, min_abund, multiplier, output):
    # open input fasta
    infile = gzip.open(input, "rt") if input.split(".")[-1] == "gz" else open(input)
    # create output fasta
    filt_outfile = open(output + ".filter.fasta", "w")
    mult_outfile = open(output + ".multiplier.fasta", "w")
    filt_mult_outfile = open(output + ".filt_mult.fasta", "w")
    # parse fasta file to extract id, multiplier, and seq
    handle = parse_fasta(infile)
    for contig_index, contig in enumerate(handle):
        contig_header, contig_id, contig_mult, contig_seq = contig
        if multiplier:
            for i in range(round(contig_mult)):
                mult_outfile.write(f">{contig_id}_{i}\n{contig_seq}\n")
        if min_abund > 0:
            if contig_mult < min_abund:
                continue
            else:
                filt_outfile.write(f">{contig_id}\n{contig_seq}\n")
                if multiplier:
                    for i in range(round(contig_mult)):
                        filt_mult_outfile.write(f">{contig_id}_{i}\n{contig_seq}\n")


def main(args=None):
    args = parse_args(args)
    multiply_sequences(args.input, float(args.abundance), args.multiplier, args.output)


if __name__ == "__main__":
    sys.exit(main())
