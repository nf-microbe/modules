#!/usr/bin/env python

import argparse
import sys

import numpy as np
import pandas as pd

from pathlib import Path


def parse_args(args=None):
    Description = "Classify sequences as complete or not based on CheckV metrics."
    Epilog = """
    Example usage:
    python completenessnfilter.py \
        -i combined_data.tsv \
        -o output.tsv
    """

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to input TSV containing combined sequence data.",
    )
    parser.add_argument(
        "-f",
        "--filters",
        help="Path to input TSV containing completeness filters to use.",
    )
    parser.add_argument(
        "-c",
        "--checkv_min_aai_completeness",
        default=80,
        help="Minimum AAI completeness to consider 'high' or 'medium' quality CheckV hits as passing.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output completeness data file in TSV format.",
    )
    return parser.parse_args(args)


def fetch_taxon(lineages: dict, lineage: str) -> str:
    if isinstance(lineage, str):
        lineage = lineage.split(";")
        for i in range(len(lineage)):
            taxon_name = ";".join(lineage[: len(lineage) - i])
            if taxon_name in lineages and lineages[taxon_name] is not None:
                return taxon_name
            else:
                return None


def main(args=None):
    args = parse_args(args)
    # load virus data
    comb_virus_data_df = pd.read_csv(args.input, sep="\t")

    if len(comb_virus_data_df) > 0:
        # load classification filters
        thresholds_df = pd.read_csv(args.filters, sep="\t")
        thresholds_dict = pd.Series(thresholds_df.min_length.values, index=thresholds_df.lineage.values).to_dict()

        # identify lowest taxonomic rank with cutoff that matches each contig
        comb_virus_data_df.loc[:, "taxon_for_threshold"] = comb_virus_data_df.apply(
            lambda row: fetch_taxon(thresholds_dict, row["lineage"]), axis=1
        )

        # map taxon information onto df
        thresholds_dict[None] = 1000000
        comb_virus_data_df["taxon_min_length"] = comb_virus_data_df.taxon_for_threshold.map(
            lambda x: thresholds_dict[x]
        )

        # classify sequence as passing or failing completeness filter
        comb_virus_data_df["completeness_filter_method"] = ""
        comb_virus_data_df["completeness_filter_method"] = np.where(
            (comb_virus_data_df["aai_confidence"].isin(["high", "medium"]))
            & (comb_virus_data_df["aai_completeness"] >= float(args.checkv_min_aai_completeness)),
            "CheckV completeness;",
            comb_virus_data_df["completeness_filter_method"],
        )
        comb_virus_data_df["completeness_filter_method"] = np.where(
            ~(comb_virus_data_df["aai_confidence"].isin(["high", "medium"]))
            & (comb_virus_data_df["contig_length"] >= comb_virus_data_df["taxon_min_length"]),
            comb_virus_data_df["completeness_filter_method"] + "taxon-specific completeness",
            comb_virus_data_df["completeness_filter_method"],
        )
        comb_virus_data_df["completeness_status"] = np.where(
            comb_virus_data_df["completeness_filter_method"] != "", "pass", "fail"
        )

        # write combined metadata to output TSV
        comb_virus_data_df.sort_values("contig_id", inplace=True)
        comb_virus_data_df.to_csv(args.output, sep="\t", index=False)

    else:
        comb_virus_data_df["completeness_filter_method"] = ""
        comb_virus_data_df["completeness_status"] = ""
        comb_virus_data_df.to_csv(args.output, sep="\t", index=False)
        comb_virus_data_df.sort_values("contig_id", inplace=True)


if __name__ == "__main__":
    sys.exit(main())
