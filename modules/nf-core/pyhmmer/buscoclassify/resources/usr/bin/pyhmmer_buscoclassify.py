#!/usr/bin/env python

import argparse
import csv
import gzip
import sys

import pandas as pd
import pyhmmer


def parse_args(args=None):
    description = "Run pyhmmer to identify extra BUSCO HMM hits."
    epilog = """
    Example usage:
    python pyhmmer.py \
        -i sequence.faa.gz \
        --archaea_hmms archaea.hmms \\
        --archaea_cutoffs archaea.cutoffs \\
        --bacteria_hmms bacteria.hmms \\
        --bacteria_cutoffs bacteria.cutoffs \\
        --output busco_hmms.tsv
        -p prefix
    """

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to input translated FASTA (gzipped or not) file containing protein sequences.",
    )
    parser.add_argument(
        "-a",
        "--archaea_hmms",
        help="Path to archaea HMM file.",
    )
    parser.add_argument(
        "-r",
        "--archaea_cutoffs",
        help="Path to archaea HMM cutoffs file.",
    )
    parser.add_argument(
        "-b",
        "--bacteria_hmms",
        help="Path to bacteria HMM file.",
    )
    parser.add_argument(
        "-c",
        "--bacteria_cutoffs",
        help="Path to bacteria HMM cutoffs file.",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        help="Output path/prefix for TSV files containing busco HMM alignment details.",
    )
    return parser.parse_args(args)


# Function to calculate hmm alignment length (coverage)
def get_hmm_coverage(domain):
    n_aligned_positions = domain.alignment.hmm_to - domain.alignment.hmm_from + 1
    return n_aligned_positions / domain.alignment.hmm_length


# Function to run pyhmmer on BUSCO HMMs
def busco_pyhmmer(faa_path, archaea_hmm, archaea_scores, bacteria_hmm, bacteria_scores, busco_hmms_out):
    # Import BUSCO cutoff scores:
    arc_cutoff_df = pd.read_csv(archaea_scores, sep="\t", index_col=False, header=None)
    arc_cutoff_dict = arc_cutoff_df.set_index(0).T.to_dict("list")

    # Import BUSCO cutoff scores:
    bact_cutoff_df = pd.read_csv(bacteria_scores, sep="\t", index_col=False, header=None)
    bact_cutoff_dict = bact_cutoff_df.set_index(0).T.to_dict("list")

    with open(busco_hmms_out, "w") as fout:
        fout.write("query\thit\tstart\tend\tevalue\tscore\thmm_coverage\n")
        with pyhmmer.easel.SequenceFile(faa_path, digital=True, alphabet=pyhmmer.easel.Alphabet.amino()) as seq_file:
            seqs = seq_file.read_block()
        with pyhmmer.plan7.HMMFile(archaea_hmm) as hmm_file:
            for hits in pyhmmer.hmmsearch(hmm_file, seqs):
                for hit in hits:
                    for domain in hit.domains.included:
                        hmm_coverage = get_hmm_coverage(domain)
                        # Store bitscore matches for each gene match. If match below cutoff and below coverage 0.5, discard.
                        if hmm_coverage >= 0.5 and hit.score >= arc_cutoff_dict[hits.query_name.decode()][0]:
                            fout.write(
                                f"{hits.query_name.decode()}\t"
                                f"{hit.name.decode()}\t"
                                f"{domain.env_from}\t"
                                f"{domain.env_to}\t"
                                f"{hit.evalue:.2E}\t"
                                f"{hit.score:.2f}\t"
                                f"{hmm_coverage:.2f}\n"
                            )
        with pyhmmer.plan7.HMMFile(bacteria_hmm) as hmm_file:
            for hits in pyhmmer.hmmsearch(hmm_file, seqs):
                for hit in hits:
                    for domain in hit.domains.included:
                        hmm_coverage = get_hmm_coverage(domain)
                        # Store bitscore matches for each gene match. If match below cutoff and below coverage 0.5, discard.
                        if hmm_coverage >= 0.5 and hit.score >= bact_cutoff_dict[hits.query_name.decode()][0]:
                            fout.write(
                                f"{hits.query_name.decode()}\t"
                                f"{hit.name.decode()}\t"
                                f"{domain.env_from}\t"
                                f"{domain.env_to}\t"
                                f"{hit.evalue:.2E}\t"
                                f"{hit.score:.2f}\t"
                                f"{hmm_coverage:.2f}\n"
                            )


def parse_hmm_hits(faa_path, busco_hmm_hits, out_path):
    data = {}
    in_faa_gunzipped = gzip.open(faa_path, "rt") if faa_path.split(".")[-1] == "gz" else open(faa_path)
    for line in in_faa_gunzipped:
        if line[0] == ">":
            contig_id = line[1:].split()[0].rpartition("_")[0]
            data[contig_id] = {
                "busco_markers": [],
            }

    busco_markers = {}
    for r in csv.DictReader(open(busco_hmm_hits), delimiter="\t"):
        busco_markers[r["hit"]] = r["query"]
    for gene, hit in busco_markers.items():
        contig_id = gene.rsplit("_", 1)[0]
        data[contig_id]["busco_markers"].append(hit)

    with open(out_path, "w") as f:
        for contig_id in data:
            n_busco = len(data[contig_id]["busco_markers"])
            busco_list = ",".join(data[contig_id]["busco_markers"])
            row = [contig_id, n_busco, busco_list]
            f.write("\t".join([str(_) for _ in row]) + "\n")


def main(args=None):
    args = parse_args(args)

    busco_pyhmmer(
        faa_path=args.input,
        archaea_hmm=args.archaea_hmms,
        archaea_scores=args.archaea_cutoffs,
        bacteria_hmm=args.bacteria_hmms,
        bacteria_scores=args.bacteria_cutoffs,
        busco_hmms_out=args.prefix + ".busco_hmms.tsv",
    )

    parse_hmm_hits(
        faa_path=args.input,
        busco_hmm_hits=args.prefix + ".busco_hmms.tsv",
        out_path=args.prefix + ".busco_markers.tsv",
    )


if __name__ == "__main__":
    sys.exit(main())
