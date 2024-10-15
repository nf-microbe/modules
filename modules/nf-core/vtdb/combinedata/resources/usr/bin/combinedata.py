#!/usr/bin/env python

import argparse
import csv
import gzip
import os
import sys

import pandas as pd
from Bio import SeqIO


def parse_args(args=None):
    description = "Combine virus data from multiple sources."
    epilog = """
    Example usage:
    python combinedata.py \
        -i sequences.fasta.gz \
        -tr trfinder_results.tsv \
        -gs genomad_aggregated_classification.tsv \
        -gg genomad_hallmarks.tsv \
        -gt genomad_taxonomy.tsv \
        -bh busco_hmm_hits.tsv \
        -ph plasmid_hmm_hits.tsv \
        -vh virus_hmm_hits.tsv \
        -cp completeness.tsv \
        -cc contamination.tsv \
        -t tantan.tsv \
        -ss nucleotides_stats.tsv \
        -o output.tsv
    """

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to input FASTA (gzipped) containing terminal repeat sequences.",
    )
    parser.add_argument(
        "-tr",
        "--trfinder",
        help="Path to input TSV file containing contigs to keep, regardless of whether they pass filters.",
    )
    parser.add_argument(
        "-gs", "--genomad_scores", help="Path to input TSV file containing geNoamd aggregated classification scores."
    )
    parser.add_argument(
        "-gg",
        "--genomad_genes",
        help="Path to input TSV file containing information about each gene detected by geNomad.",
    )
    parser.add_argument(
        "-gt",
        "--genomad_taxa",
        help="Path to input TSV file containing geNomad's marker-gene taxonomy.",
    )
    parser.add_argument(
        "-bh",
        "--busco_hmms",
        help="Path to input TSV file containing BUSCO HMM hits.",
    )
    parser.add_argument(
        "-ph",
        "--plasmid_hmms",
        help="Path to input TSV file containing plasmid HMM hits.",
    )
    parser.add_argument(
        "-vh",
        "--virus_hmms",
        help="Path to input TSV file containing virus HMM hits.",
    )
    parser.add_argument(
        "-cp",
        "--completeness",
        help="Path to input TSV file containing CheckV's completeness output.",
    )
    parser.add_argument(
        "-cc",
        "--contamination",
        help="Path to input TSV file containing CheckV's contamination output.",
    )
    parser.add_argument(
        "-t",
        "--tantan",
        help="Path to input TSV file containing tantan's output.",
    )
    parser.add_argument(
        "-ss",
        "--sequence_stats",
        help="Path to input TSV file containing nucleotide stats (GC %, CDS density, and 'N' percentage).",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Location to output combined CSV file to.",
    )
    return parser.parse_args(args)


def combine_virus_data(
    input_fasta,
    output,
    trfinder=None,
    genomad_scores=None,
    genomad_genes=None,
    genomad_taxa=None,
    busco_hmms=None,
    plasmid_hmms=None,
    virus_hmms=None,
    completeness=None,
    contamination=None,
    tantan=None,
    sequence_stats=None,
):
    # load trfinder data
    if trfinder:
        if os.path.getsize(trfinder) > 0:
            trfinder_df = pd.read_csv(
                trfinder,
                sep="\t",
                header=0,
                index_col="contig_name",
                usecols=[
                    "contig_name",
                    "contig_len",
                    "tr_type",
                    "tr_seq",
                    "tr_len",
                    "tr_nt_acgt_count",
                    "tr_nt_n_count",
                    "tr_nt_max_freq",
                ],
            )
    if "trfinder_df" not in locals():
        trfinder_df = pd.DataFrame(
            columns=["contig_len", "tr_type", "tr_seq", "tr_len", "tr_nt_acgt_count", "tr_nt_n_count", "tr_nt_max_freq"]
        )

    # load genomad agg scores
    if genomad_scores:
        if os.path.getsize(genomad_scores) > 0:
            genomad_scores_df = pd.read_csv(
                genomad_scores,
                sep="\t",
                header=0,
                index_col="seq_name",
            )
    if "genomad_scores_df" not in locals():
        genomad_scores_df = pd.DataFrame(columns=["chromosome_score", "plasmid_score", "virus_score"])

    # load genomad marker gene data
    if genomad_genes:
        if os.path.getsize(genomad_genes) > 0:
            # convert genomad genes to marker stats
            genomad_marker_stats = {}
            input_gunzipped = gzip.open(input_fasta, "rt") if input_fasta.split(".")[-1] == "gz" else open(input_fasta)
            for r in SeqIO.parse(input_gunzipped, "fasta"):
                genomad_marker_stats[r.id] = {"contig_id": r.id, "uscg_count": 0, "plasmid_count": 0, "virus_count": 0}
            # sum counts across genes within contigs
            for r in csv.DictReader(open(genomad_genes), delimiter="\t"):
                contig_id = r["gene"].rsplit("_", 1)[0]
                genomad_marker_stats[contig_id]["uscg_count"] += int(r["uscg"])
                genomad_marker_stats[contig_id]["plasmid_count"] += int(r["plasmid_hallmark"])
                genomad_marker_stats[contig_id]["virus_count"] += int(r["virus_hallmark"])
            # convert dict to df
            if genomad_marker_stats:
                genomad_hallmarks_df = pd.DataFrame.from_dict(genomad_marker_stats, orient="index").drop(
                    "contig_id", axis=1
                )
    if "genomad_hallmarks_df" not in locals():
        genomad_hallmarks_df = pd.DataFrame(columns=["uscg_count", "plasmid_count", "virus_count"])

    # load genoamd taxonomy data
    if genomad_taxa:
        if os.path.getsize(genomad_taxa) > 0:
            genomad_taxa_df = pd.read_csv(
                genomad_taxa, sep="\t", header=0, index_col="seq_name", usecols=["seq_name", "lineage"]
            )
    if "genomad_taxa_df" not in locals():
        genomad_taxa_df = pd.DataFrame(columns=["lineage"])

    # load busco HMM data
    if busco_hmms:
        if os.path.getsize(busco_hmms) > 0:
            busco_hmms_df = pd.read_csv(
                busco_hmms,
                sep="\t",
                header=None,
                names=["contig_id", "busco_count", "busco_list"],
                index_col="contig_id",
                usecols=["contig_id", "busco_count", "busco_list"],
            )
    if "busco_hmms_df" not in locals():
        busco_hmms_df = pd.DataFrame(columns=["busco_count", "busco_list"])

    # load plasmid HMM data
    if plasmid_hmms:
        if os.path.getsize(plasmid_hmms) > 0:
            plasmid_hmms_df = pd.read_csv(
                plasmid_hmms,
                sep="\t",
                header=None,
                names=["contig_id", "plasmid_count", "plasmid_list"],
                index_col="contig_id",
                usecols=["contig_id", "plasmid_count", "plasmid_list"],
            )
    if "plasmid_hmms_df" not in locals():
        plasmid_hmms_df = pd.DataFrame(columns=["plasmid_count", "plasmid_list"])

    # load busco HMM data
    if virus_hmms:
        if os.path.getsize(virus_hmms) > 0:
            virus_hmms_df = pd.read_csv(
                virus_hmms,
                sep="\t",
                header=None,
                names=["contig_id", "virus_count", "virus_list"],
                index_col="contig_id",
                usecols=["contig_id", "virus_count", "virus_list"],
            )
    if "virus_hmms_df" not in locals():
        virus_hmms_df = pd.DataFrame(columns=["virus_count", "virus_list"])

    # load checkv completeness data
    if completeness:
        if os.path.getsize(completeness) > 0:
            completeness_df = pd.read_csv(
                completeness,
                sep="\t",
                header=0,
                index_col="contig_id",
                usecols=[
                    "contig_id",
                    "viral_length",
                    "aai_completeness",
                    "aai_confidence",
                    "hmm_completeness_lower",
                    "hmm_num_hits",
                    "kmer_freq",
                ],
            )
            completeness_df.hmm_completeness_lower = completeness_df.hmm_completeness_lower.round(4)
    if "completeness_df" not in locals():
        completeness_df = pd.DataFrame(
            columns=["aai_completeness", "aai_confidence", "hmm_completeness_lower", "hmm_num_hits", "kmer_freq"]
        )

    # load checkv contamination data
    if contamination:
        if os.path.getsize(contamination) > 0:
            contamination_df = pd.read_csv(
                contamination,
                sep="\t",
                header=0,
                index_col="contig_id",
                usecols=["contig_id", "viral_genes", "host_genes"],
            )
    if "contamination_df" not in locals():
        contamination_df = pd.DataFrame(columns=["viral_genes", "host_genes"])

    # load tantan data
    if tantan:
        if os.path.getsize(tantan) > 0:
            input_gunzipped = gzip.open(input_fasta, "rt") if input_fasta.split(".")[-1] == "gz" else open(input_fasta)
            tantan_dict = {}
            for r in SeqIO.parse(input_gunzipped, "fasta"):
                tantan_dict[r.id] = 0
            # tantan_dict = dict([[contig_id, 0] for contig_id in genomad_scores_df.index])
            with open(tantan) as tantan_file:
                for line in tantan_file:
                    contig_id, start, end = line.split()
                    length = int(end) - int(start) + 1
                    tantan_dict[contig_id] += length
            if tantan_dict:
                tantan_df = pd.DataFrame.from_dict(tantan_dict, orient="index")
                tantan_df.columns = ["tantan_len"]
    if "tantan_df" not in locals():
        tantan_df = pd.DataFrame(columns=["tantan_len"])

    # load nucleotide stats
    if sequence_stats:
        sequence_stats_df = pd.read_csv(
            sequence_stats,
            sep="\t",
            header=0,
            index_col="contig_id",
        )
    if "sequence_stats_df" not in locals():
        sequence_stats_df = pd.DataFrame(
            columns=["contig_length", "cds_length", "cds_density", "gc_percent", "n_count"]
        )

    # combine quality data from all sources by contig_id (index)
    comb_virus_data_df = pd.concat(
        [
            trfinder_df,
            genomad_scores_df,
            genomad_hallmarks_df,
            genomad_taxa_df,
            busco_hmms_df,
            plasmid_hmms_df,
            virus_hmms_df,
            completeness_df,
            contamination_df,
            tantan_df,
            sequence_stats_df,
        ],
        axis=1,
    )

    # perform calculations
    if tantan and sequence_stats:
        comb_virus_data_df.insert(
            29, "tantan_freq", comb_virus_data_df["tantan_len"] / comb_virus_data_df["contig_length"]
        )
        comb_virus_data_df.tantan_freq = comb_virus_data_df.tantan_freq.astype(float)
        comb_virus_data_df.tantan_freq = comb_virus_data_df.tantan_freq.round(4)
    if sequence_stats:
        comb_virus_data_df.insert(
            len(comb_virus_data_df.columns),
            "n_freq",
            comb_virus_data_df["n_count"] / comb_virus_data_df["contig_length"],
        )

    # write combined metadata to output TSV
    comb_virus_data_df.insert(0, "contig_id", comb_virus_data_df.index)
    comb_virus_data_df['contig_id'] = comb_virus_data_df['contig_id'].astype(str)
    comb_virus_data_df.sort_values("contig_id", inplace=True)
    comb_virus_data_df.to_csv(output, sep="\t", index=False)


def main(args=None):
    args = parse_args(args)
    combine_virus_data(
        args.input,
        args.output,
        args.trfinder,
        args.genomad_scores,
        args.genomad_genes,
        args.genomad_taxa,
        args.busco_hmms,
        args.plasmid_hmms,
        args.virus_hmms,
        args.completeness,
        args.contamination,
        args.tantan,
        args.sequence_stats,
    )


if __name__ == "__main__":
    sys.exit(main())
