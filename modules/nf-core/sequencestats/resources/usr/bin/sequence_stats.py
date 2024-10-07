#!/usr/bin/env python

import argparse
import gzip
import sys

import Bio.SeqIO


def parse_args(args=None):
    description = "Calculate nucleotide stats (GC %, CDS density, N count) on input FASTA."
    epilog = "Example usage: python sequence_stats.py -i sequences.fasta.gz -p proteins.faa.gz -o nucleotide_stat.tsv"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to input nucleotide fasta to be used in calculating sequence statistics.",
    )
    parser.add_argument(
        "-p",
        "--proteins",
        help="Path to input FAA file (gzipped) containing protein sequences associated with input FASTA.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output TSV file containing nucleotide statistics.",
    )
    return parser.parse_args(args)


def calc_gc(x):
    x = x.upper()
    gc = round(100.0 * (x.count("G") + x.count("C")) / len(x), 2)
    return gc


def nucleotide_stats(input, proteins, output):
    stats = {}
    input_gunzipped = gzip.open(input, "rt") if input.split(".")[-1] == "gz" else open(input)
    for r in Bio.SeqIO.parse(input_gunzipped, "fasta"):
        stats[r.id] = {
            "contig_id": r.id,
            "contig_length": len(r.seq),
            "cds_length": 0,
            "gc_content": calc_gc(str(r.seq)),
            "n_count": str(r.seq).upper().count("N"),
        }

    proteins_gunzipped = gzip.open(proteins, "rt") if input.split(".")[-1] == "gz" else open(proteins)
    for line in proteins_gunzipped:
        if line[0] == ">":
            r = line.split()
            contig_id = r[0][1:].rsplit("_", 1)[0]
            stats[contig_id]["cds_length"] += int(r[4]) - int(r[2]) + 1

    with open(output, "w") as out:
        fields = ["contig_id", "contig_length", "cds_length", "cds_density", "gc_content", "n_count"]
        out.write("\t".join(fields) + "\n")
        for r in stats.values():
            r["cds_density"] = round(100 * r["cds_length"] / r["contig_length"], 2)
            row = [str(r[f]) for f in fields]
            out.write("\t".join(row) + "\n")


def main(args=None):
    args = parse_args(args)
    nucleotide_stats(args.input, args.proteins, args.output)


if __name__ == "__main__":
    sys.exit(main())
