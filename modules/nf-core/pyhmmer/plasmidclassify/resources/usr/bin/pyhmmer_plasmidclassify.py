#!/usr/bin/env python

import argparse
import csv
import gzip
import sys

import pyhmmer


def parse_args(args=None):
    description = "Run pyhmmer to identify extra plasmid HMM hits."
    epilog = """
    Example usage:
    python pyhmmer_plasmidclassify.py \
        -i proteins.faa.gz \
        -l plasmid_hallmarks.hmm \
        -p prefix
    """

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to input translated FASTA (gzipped or not) file containing protein sequences.",
    )
    parser.add_argument(
        "-l",
        "--plasmid_hallmarks",
        help="Path to plasmid hallmarks hmm file.",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        help="Output path/prefix for TSV files containing plasmid HMM alignment details.",
    )
    return parser.parse_args(args)


# Function to calculate hmm alignment length (coverage)
def get_hmm_coverage(domain):
    n_aligned_positions = domain.alignment.hmm_to - domain.alignment.hmm_from + 1
    return n_aligned_positions / domain.alignment.hmm_length


# Function to run pyhmmer on plasmid HMMs
def plasmid_pyhmmer(faa_path, plasmid_hallmarks_hmm, plasmid_hmms_out):
    with open(plasmid_hmms_out, "w") as fout:
        fout.write("query\thit\tstart\tend\tevalue\tscore\thmm_coverage\n")
        with pyhmmer.easel.SequenceFile(faa_path, digital=True, alphabet=pyhmmer.easel.Alphabet.amino()) as seq_file:
            seqs = seq_file.read_block()
        with pyhmmer.plan7.HMMFile(plasmid_hallmarks_hmm) as hmm_file:
            for hits in pyhmmer.hmmsearch(hmm_file, seqs, bit_cutoffs="gathering"):
                for hit in hits:
                    for domain in hit.domains.included:
                        hmm_coverage = get_hmm_coverage(domain)
                        if hmm_coverage >= 0.5:
                            fout.write(
                                f"{hits.query_name.decode()}\t"
                                f"{hit.name.decode()}\t"
                                f"{domain.env_from}\t"
                                f"{domain.env_to}\t"
                                f"{hit.evalue:.2E}\t"
                                f"{hit.score:.2f}\t"
                                f"{hmm_coverage:.2f}\n"
                            )


def parse_hmm_hits(faa_path, plasmid_hmm_hits, out_path):
    data = {}
    in_faa_gunzipped = gzip.open(faa_path, "rt") if faa_path.split(".")[-1] == "gz" else open(faa_path)
    for line in in_faa_gunzipped:
        if line[0] == ">":
            contig_id = line[1:].split()[0].rpartition("_")[0]
            data[contig_id] = {
                "plasmid_markers": [],
            }

    plasmid_markers = {}
    for r in csv.DictReader(open(plasmid_hmm_hits), delimiter="\t"):
        plasmid_markers[r["hit"]] = r["query"]
    for gene, hit in plasmid_markers.items():
        contig_id = gene.rsplit("_", 1)[0]
        data[contig_id]["plasmid_markers"].append(hit)

    with open(out_path, "w") as f:
        for contig_id in data:
            n_plasmid = len(data[contig_id]["plasmid_markers"])
            plasmid_list = ",".join(data[contig_id]["plasmid_markers"])
            row = [contig_id, n_plasmid, plasmid_list]
            f.write("\t".join([str(_) for _ in row]) + "\n")


def main(args=None):
    args = parse_args(args)

    plasmid_pyhmmer(
        faa_path=args.input,
        plasmid_hallmarks_hmm=args.plasmid_hallmarks,
        plasmid_hmms_out=args.prefix + ".plasmid_hmms.tsv",
    )

    parse_hmm_hits(
        faa_path=args.input,
        plasmid_hmm_hits=args.prefix + ".plasmid_hmms.tsv",
        out_path=args.prefix + ".plasmid_markers.tsv",
    )


if __name__ == "__main__":
    sys.exit(main())
