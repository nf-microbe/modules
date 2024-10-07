#!/usr/bin/env python

import argparse
import csv
import gzip
import sys

import pyhmmer


def parse_args(args=None):
    description = "Run pyhmmer to identify extra viral HMM hits."
    epilog = """
    Example usage:
    python pyhmmervirus.py \
        -i proteins.faa.gz \
        -d DJR_MCP_virus_hallmarks.hmm \
        -n Inovirus_MCP_virus_hallmarks.hmm \
        -l pleolipoviridae.hmm \
        -p prefix
    """

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to input translated FASTA (gzipped or not) file containing protein sequences.",
    )
    parser.add_argument(
        "-d",
        "--djr",
        help="Path to DJR MCP hmm file.",
    )
    parser.add_argument(
        "-n",
        "--inovirus",
        help="Path to Inovirus MCP hmm file.",
    )
    parser.add_argument(
        "-l",
        "--pleolipoviridae",
        help="Path to pleolipoviridae hmm file.",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        help="Output path/prefix for TSV files containing virus HMM alignment details.",
    )
    return parser.parse_args(args)


# Function to calculate hmm alignment length (coverage)
def get_hmm_coverage(domain):
    n_aligned_positions = domain.alignment.hmm_to - domain.alignment.hmm_from + 1
    return n_aligned_positions / domain.alignment.hmm_length


# Function to run pyhmmer on viral HMMs
def virus_pyhmmer(faa_path, djr_hmm, inovirus_hmm, pleolipoviridae_hmm, virus_hmms_out):
    with open(virus_hmms_out, "w") as fout:
        fout.write("query\thit\tstart\tend\tevalue\tscore\thmm_coverage\n")
        with pyhmmer.easel.SequenceFile(faa_path, digital=True, alphabet=pyhmmer.easel.Alphabet.amino()) as seq_file:
            seqs = seq_file.read_block()
        with pyhmmer.plan7.HMMFile(djr_hmm) as hmm_file:
            for hits in pyhmmer.hmmsearch(hmm_file, seqs, Z=10_000, E=1e-12):
                for hit in hits:
                    for domain in hit.domains:
                        if not domain.included:
                            continue
                        if (hmm_coverage := get_hmm_coverage(domain)) >= 0.7:
                            fout.write(
                                f"{hits.query_name.decode()}\t"
                                f"{hit.name.decode()}\t"
                                f"{domain.env_from}\t"
                                f"{domain.env_to}\t"
                                f"{hit.evalue:.2E}\t"
                                f"{hit.score:.2f}\t"
                                f"{hmm_coverage:.2f}\n"
                            )
        with pyhmmer.plan7.HMMFile(inovirus_hmm) as hmm_file:
            for hits in pyhmmer.hmmsearch(hmm_file, seqs, Z=10_000, E=6e-05):
                for hit in hits:
                    for domain in hit.domains:
                        if not domain.included:
                            continue
                        if (hmm_coverage := get_hmm_coverage(domain)) >= 0.8:
                            fout.write(
                                f"{hits.query_name.decode()}\t"
                                f"{hit.name.decode()}\t"
                                f"{domain.env_from}\t"
                                f"{domain.env_to}\t"
                                f"{hit.evalue:.2E}\t"
                                f"{hit.score:.2f}\t"
                                f"{hmm_coverage:.2f}\n"
                            )
        with pyhmmer.plan7.HMMFile(pleolipoviridae_hmm) as hmm_file:
            for hits in pyhmmer.hmmsearch(hmm_file, seqs, Z=10_000, E=1e-9):
                for hit in hits:
                    for domain in hit.domains:
                        if not domain.included:
                            continue
                        if (hmm_coverage := get_hmm_coverage(domain)) >= 0.5:
                            fout.write(
                                f"{hits.query_name.decode()}\t"
                                f"{hit.name.decode()}\t"
                                f"{domain.env_from}\t"
                                f"{domain.env_to}\t"
                                f"{hit.evalue:.2E}\t"
                                f"{hit.score:.2f}\t"
                                f"{hmm_coverage:.2f}\n"
                            )


def parse_hmm_hits(faa_path, virus_hmm_hits, out_path):
    data = {}
    in_faa_gunzipped = gzip.open(faa_path, "rt") if faa_path.split(".")[-1] == "gz" else open(faa_path)
    for line in in_faa_gunzipped:
        if line[0] == ">":
            contig_id = line[1:].split()[0].rpartition("_")[0]
            data[contig_id] = {
                "virus_markers": [],
            }

    virus_markers = {}
    for r in csv.DictReader(open(virus_hmm_hits), delimiter="\t"):
        virus_markers[r["hit"]] = r["query"]
    for gene, hit in virus_markers.items():
        contig_id = gene.rsplit("_", 1)[0]
        data[contig_id]["virus_markers"].append(hit)

    with open(out_path, "w") as f:
        for contig_id in data:
            n_virus = len(data[contig_id]["virus_markers"])
            virus_list = ",".join(data[contig_id]["virus_markers"])
            row = [contig_id, n_virus, virus_list]
            f.write("\t".join([str(_) for _ in row]) + "\n")


def main(args=None):
    args = parse_args(args)

    virus_pyhmmer(
        faa_path=args.input,
        djr_hmm=args.djr,
        inovirus_hmm=args.inovirus,
        pleolipoviridae_hmm=args.pleolipoviridae,
        virus_hmms_out=args.prefix + ".virus_hmms.tsv",
    )

    parse_hmm_hits(
        faa_path=args.input,
        virus_hmm_hits=args.prefix + ".virus_hmms.tsv",
        out_path=args.prefix + ".virus_markers.tsv",
    )


if __name__ == "__main__":
    sys.exit(main())
