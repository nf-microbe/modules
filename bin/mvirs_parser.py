#!/usr/bin/env python

import argparse
import gzip

import Bio.SeqIO
import pandas as pd


def find_direct_repeat(contig_seq, mge_start, mge_end, max_repeat=30, max_att=100):
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

    longest_dr = ""
    mge_new_start = mge_start
    mge_new_end = mge_end

    # Start at mvirs edges, and step inwards (towards MGE center) to look for repeats
    for start_shift in range(0, max_repeat):
        # Extract the flanks
        left_flank = contig_seq[: mge_start + start_shift]  # Left side of the subsequence
        right_flank = contig_seq[mge_end - start_shift :]  # Right side of the subsequence

        # Start with the largest possible repeat size and work down
        for dr_len in range(max_repeat, 0, -1):
            # Make sure enough flanking sequence is present for the DR + att site
            if dr_len + max_att > len(left_flank) or dr_len + max_att > len(right_flank):
                continue

            # Get the potential repeat on both flanks
            left_repeat = left_flank[-dr_len:]
            right_repeat = right_flank[:dr_len]

            # Check if the repeat is the same on both sides
            if left_repeat == right_repeat:
                # Make sure repeat is adjacent to mvirs edge
                if dr_len >= start_shift:
                    # Check if repeat is longer that previous best
                    if dr_len > len(longest_dr):
                        mge_new_start = mge_start + start_shift
                        mge_new_end = mge_end - start_shift
                        longest_dr = left_repeat

    return mge_new_start, mge_new_end, len(longest_dr), longest_dr


def find_att_sites(contig_seq, mge_start, mge_end, dr_seq, att_len, context_len):
    b1 = contig_seq[: mge_start - len(dr_seq)][-att_len:]
    b2 = contig_seq[mge_end + len(dr_seq) :][:att_len]
    attb = b1 + dr_seq + b2
    con1 = contig_seq[: mge_start - len(dr_seq)][-(context_len + att_len) :]
    con2 = contig_seq[mge_end + len(dr_seq) :][: (context_len + att_len)]
    context = con1 + dr_seq + con2
    p1 = contig_seq[mge_start : mge_start + att_len]
    p2 = contig_seq[mge_end - att_len : mge_end]
    attp = p2 + dr_seq + p1

    return attb, context, attp


def make_integrase_output(faa_path, domtbl_path):
    rows = []
    proteins_open = gzip.open(faa_path, "rt") if faa_path.split(".")[-1] == "gz" else open(faa_path)
    proteins = dict([[r.id, r.seq] for r in Bio.SeqIO.parse(proteins_open, "fasta")])
    domtbl_path_open = gzip.open(domtbl_path, "rt") if domtbl_path.split(".")[-1] == "gz" else open(domtbl_path)
    for line in domtbl_path_open:
        if line[0] == "#":
            continue
        parts = line.split()
        start, stop, strand_flag, gene_id, hmm_id, evalue, score = (
            parts[-8:][1],
            parts[-8:][3],
            parts[-8:][5],
            parts[0],
            parts[3],
            parts[4],
            parts[5],
        )
        strand = "+" if strand_flag == "1" else "-"
        mge_id, cds_num = gene_id.rsplit("_", 1)
        protein = proteins[gene_id]
        row = [mge_id, cds_num, hmm_id, start, stop, strand, score, evalue, protein]
        rows.append(row)
    df = pd.DataFrame(rows)
    if len(df.columns) == 9:
        df.columns = ["contig_id", "cds_num", "hmm_id", "start", "stop", "strand", "score", "evalue", "protein"]
    else:
        df.columns = ["contig_id", "cds_num", "hmm_id", "start", "stop", "strand", "score", "evalue", "protein"]

    return df


def find_direct_repeats(fna_path, mvirs_path, faa_path, domtbl_path, summary_path, max_repeat, att_len, context_len):
    """
    Process MGEs sequences to find direct repeats flanking mobile genetic elements (MGEs).

    Args:
        fna_path (str): Path to the compressed contigs file (.fna.gz).
        mvirs_path (str): Path to the compressed prophages file (.fasta.gz).
        out_path (str): Path to save the output TSV file containing direct repeat information.
        max_repeat (int): Maximum repeat size to search for (default is 30).

    The function parses the DNA sequences from the provided input files, extracts the
    regions of interest, searches for the largest direct repeat flanking each prophage,
    and writes the results to a TSV file.
    """

    prophages_open = gzip.open(mvirs_path, "rt") if mvirs_path.split(".")[-1] == "gz" else open(mvirs_path)
    prophages = dict([[r.id, r.seq] for r in Bio.SeqIO.parse(prophages_open, "fasta")])
    contigs_open = gzip.open(fna_path, "rt") if fna_path.split(".")[-1] == "gz" else open(fna_path)
    contigs = dict([[r.id, r.seq] for r in Bio.SeqIO.parse(contigs_open, "fasta")])

    rows = []
    for mge_id, prophage_seq in prophages.items():
        assembly_id = mge_id.split("_")[0].split(":")[0]
        contig_id = mge_id.split(":")[0]
        mge_start, mge_end = (int(_) for _ in mge_id.split(":")[-1].split("-"))
        opr = int(mge_id.split(" ")[1].split("=")[1].split('-')[0])
        sr = int(mge_id.split(" ")[1].split("=")[2].split('-')[0])
        contig_seq = contigs[contig_id]
        mge_new_start, mge_new_end, dr_len, dr_seq = find_direct_repeat(contig_seq, mge_start, mge_end, max_repeat)
        attb, context, attp = find_att_sites(contig_seq, mge_new_start, mge_new_end, dr_seq, att_len, context_len)
        row = [mge_id, assembly_id, contig_id, mge_new_start, mge_new_end, dr_len, opr, sr, dr_seq, attb, context, attp]
        rows.append(row)

    dr_df = pd.DataFrame(rows)
    dr_df.columns = [
        "mge_id",
        "assembly_id",
        "contig_id",
        "mge_start",
        "mge_end",
        "opr_count",
        "sr_count",
        "dr_len",
        "dr_seq",
        "attb",
        "context",
        "attp",
    ]
    integrase_df = make_integrase_output(faa_path, domtbl_path)
    summary_df = dr_df.merge(integrase_df, on="contig_id", how="outer")
    summary_df.to_csv(summary_path, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse mVIRs output to identify recombinases, find direct repeats, and att sites."
    )
    parser.add_argument("--mvirs", type=str, help="Path to the prophages .fasta.gz file")
    parser.add_argument("--fna", type=str, help="Path to the contigs .fna.gz file")
    parser.add_argument("--faa", type=str, help="Path to the protein .faa.gz file")
    parser.add_argument("--hmmsearch", type=str, help="Path to the recombinase HMMsearch output file")
    parser.add_argument(
        "--att_len", metavar="INT", type=int, default=100, help="Length for attachment sites (default: 100)"
    )
    parser.add_argument(
        "--context_len", metavar="INT", type=int, default=2500, help="Length for attachment sites (default: 2500)"
    )
    parser.add_argument(
        "--max_repeat", metavar="INT", type=int, default=30, help="Maximum repeat size to search for (default: 50)"
    )
    parser.add_argument("--prefix", type=str, help="Prefix for output files")

    args = parser.parse_args()

    args.summary_path = f"{args.prefix}.summary.tsv"

    print("Finding direct repeats and making summary file")
    find_direct_repeats(
        args.fna,
        args.mvirs,
        args.faa,
        args.hmmsearch,
        args.summary_path,
        args.max_repeat,
        args.att_len,
        args.context_len,
    )
