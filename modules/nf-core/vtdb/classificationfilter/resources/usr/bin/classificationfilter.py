#!/usr/bin/env python

import argparse
import sys

import numpy as np
import pandas as pd

from pathlib import Path


def parse_args(args=None):
    Description = "Classify sequences as viral/non-viral based on tool metrics."
    Epilog = """
    Example usage:
    python classificationfilter.py \
        -i combined_data.tsv \
        -f classification_filters.tsv \
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
        help="Path to input TSV containing classification filters to use.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output classification data file in TSV format.",
    )
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)
    # load virus data
    comb_virus_data_df = pd.read_csv(args.input, sep="\t")

    if len(comb_virus_data_df) > 0:
        # load classification filters
        thresholds_df = pd.read_csv(args.filters, sep="\t")
        viral_thresholds_df = thresholds_df[thresholds_df["classification"] == "viral"]
        nonviral_thresholds_df = thresholds_df[thresholds_df["classification"] == "non-viral"]

        # iterate through viral filters
        comb_virus_data_df["viral_classification_method"] = ""
        for index, row in viral_thresholds_df.iterrows():
            print(row["criteria"])
            comb_virus_data_df.iloc[
                comb_virus_data_df.query(row["criteria"]).index,
                len(comb_virus_data_df.columns) - 1,
            ] = (
                comb_virus_data_df.iloc[
                    comb_virus_data_df.query(row["criteria"]).index,
                    len(comb_virus_data_df.columns) - 1,
                ]
                + row["name"]
                + ";"
            )

        # iterate through non-viral filters
        comb_virus_data_df["nonviral_classification_method"] = ""
        for index, row in nonviral_thresholds_df.iterrows():
            comb_virus_data_df.iloc[
                comb_virus_data_df.query(row["criteria"]).index,
                len(comb_virus_data_df.columns) - 1,
            ] = (
                comb_virus_data_df.iloc[
                    comb_virus_data_df.query(row["criteria"]).index,
                    len(comb_virus_data_df.columns) - 1,
                ]
                + row["name"]
                + ";"
            )

        # classify sequence as viral or non-viral
        comb_virus_data_df["virus_classification"] = "non-viral"
        comb_virus_data_df["virus_classification"] = np.where(
            (comb_virus_data_df["viral_classification_method"] != "")
            & (comb_virus_data_df["nonviral_classification_method"] == ""),
            "viral",
            "non-viral",
        )

        # write combined metadata to output TSV
        comb_virus_data_df.sort_values("contig_id", inplace=True)
        comb_virus_data_df.to_csv(args.output, sep="\t", index=False)

    else:
        comb_virus_data_df["viral_classification_method"] = ""
        comb_virus_data_df["virus_classification"] = ""
        comb_virus_data_df.sort_values("contig_id", inplace=True)
        comb_virus_data_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    sys.exit(main())
