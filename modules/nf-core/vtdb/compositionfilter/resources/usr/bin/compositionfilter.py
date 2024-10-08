#!/usr/bin/env python

import argparse
import sys

import numpy as np
import pandas as pd

from pathlib import Path


def parse_args(args=None):
    Description = "Filter sequences based on previously calculated composition metrics."
    Epilog = """
    Example usage:
    python compositionfilter.py \
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
        "-l",
        "--min_len",
        help="Minimum sequence length (used as lower threshold for first quantile).",
        default=1500,
        type=int,
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


def filter_by_composition(comb_data: pd.DataFrame, min_len: int, max_len: int, filters: pd.DataFrame) -> pd.DataFrame:
    print(min_len, max_len)
    len_filt_df = comb_data.query(f"`contig_length` >= {min_len} and `contig_length` < {max_len}").reset_index(
        drop=True
    )
    if len(len_filt_df) > 0:
        for index, row in filters.query(f"`genome_length` == {max_len}").iterrows():
            len_filt_df.iloc[len_filt_df.query(row["criteria"]).index, len(comb_data.columns) - 1] = (
                len_filt_df.iloc[len_filt_df.query(row["criteria"]).index, len(comb_data.columns) - 1]
                + row["criteria"]
                + ";"
            )
        return len_filt_df


def main(args=None):
    args = parse_args(args)
    # load virus data
    comb_virus_data_df = pd.read_csv(args.input, sep="\t")

    if len(comb_virus_data_df) > 0:
        # load compoisition filters
        thresholds_df = pd.read_csv(args.filters, sep="\t")

        # run composition filters on each length quantile and concatenate
        genome_length_thresholds = [args.min_len] + thresholds_df["genome_length"].unique().tolist()
        comb_virus_data_df["composition_filter_method"] = ""
        comb_compos_df = pd.concat(
            [
                filter_by_composition(
                    comb_virus_data_df, genome_length_thresholds[i], genome_length_thresholds[i + 1], thresholds_df
                )
                for i in range(len(genome_length_thresholds) - 1)
            ]
        )

        # classify sequence as passing or failing
        comb_compos_df["composition_status"] = "passing"
        comb_compos_df["composition_status"] = np.where(
            comb_compos_df["composition_filter_method"] == "", "pass", "fail"
        )

        # write combined metadata to output TSV
        comb_virus_data_df.sort_values("contig_id", inplace=True)
        comb_compos_df.to_csv(args.output, sep="\t", index=False)

    else:
        comb_virus_data_df["composition_filter_method"] = ""
        comb_virus_data_df["composition_status"] = ""
        comb_virus_data_df.sort_values("contig_id", inplace=True)
        comb_virus_data_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    sys.exit(main())
