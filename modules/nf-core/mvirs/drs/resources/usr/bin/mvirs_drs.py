#!/usr/bin/env python

import argparse
import gzip

import Bio.SeqIO
import pandas as pd


def find_direct_repeat(seq, start, end, max_repeat=30):
    """
    Find the maximum direct repeat flanking a subsequence within a larger DNA sequence.

    Args:
        seq (str): The larger DNA sequence.
        start (int): Start index of the subseq.
        end (int): Start index of the subseq.

    Returns:
        int: The length of the maximum direct repeat.
        str: The repeat sequence itself.
    """

    # Extract the flanks
    left_flank = seq[:start]  # Left side of the subsequence
    right_flank = seq[end:]  # Right side of the subsequence

    # Start with the largest possible repeat size and work down
    for repeat_size in range(max_repeat, 0, -1):
        # Make sure enough flanking sequence is present
        if repeat_size > len(left_flank) or repeat_size > len(right_flank):
            continue

        # Get the potential repeat on both flanks
        left_repeat = left_flank[-repeat_size:]
        right_repeat = right_flank[:repeat_size]

        # Check if the repeat is the same on both sides
        if left_repeat == right_repeat:
            return repeat_size, left_repeat

    # If no repeat is found, return 0 and an empty string
    return 0, ""


def main(fna_path, mvirs_path, out_path, max_repeat):
    """
    Process DNA sequences to find direct repeats flanking mobile genetic elements (MGEs).

    Args:
        fna_path (str): Path to the compressed contigs file (.fna.gz).
        mvirs_path (str): Path to the compressed prophages file (.fasta.gz).
        out_path (str): Path to save the output TSV file containing direct repeat information.
        max_repeat (int): Maximum repeat size to search for (default is 30).
        reformat_ids (bool): Whether to reformat contig IDs by removing leading zeros.

    The function parses the DNA sequences from the provided input files, extracts the
    regions of interest, searches for the largest direct repeat flanking each prophage,
    and writes the results to a TSV file.
    """

    prophages = dict([[r.id, r.seq] for r in Bio.SeqIO.parse(gzip.open(mvirs_path, "rt"), "fasta")])
    contigs = dict([[r.id, r.seq] for r in Bio.SeqIO.parse(gzip.open(fna_path, "rt"), "fasta")])

    rows = []
    for prophage_id, prophage_seq in prophages.items():
        assembly_id = prophage_id.split("_")[0].split(":")[0]
        contig_id = prophage_id.split(":")[0]
        mge_start, mge_end = [int(_) for _ in prophage_id.split(":")[-1].split("-")]
        contig_seq = contigs[contig_id]
        dr_len, dr_seq = find_direct_repeat(contig_seq, mge_start, mge_end, max_repeat)
        row = [prophage_id, assembly_id, contig_id, mge_start, mge_end, dr_len, dr_seq]
        rows.append(row)

    df = pd.DataFrame(rows)
    df.columns = ["prophage_id", "assembly_id", "contig_id", "mge_start", "mge_end", "dr_len", "dr_seq"]
    df.to_csv(out_path, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find direct repeats flanking prophages in a DNA sequence.")
    parser.add_argument("fna_path", type=str, help="Path to the contigs .fna.gz file")
    parser.add_argument("mvirs_path", type=str, help="Path to the prophages .fasta.gz file")
    parser.add_argument("out_path", type=str, help="Output file path to save direct repeat information")
    parser.add_argument(
        "--max_repeat", metavar="INT", type=int, default=30, help="Maximum repeat size to search for (default: 30)"
    )

    args = parser.parse_args()

    main(args.fna_path, args.mvirs_path, args.out_path, args.max_repeat)
