#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Erin Nicole Decocker
# edecocke@charlotte.edu
# ID: 801442694
# AI usage acknowledgment in READme file.

"""
main.py

This script serves as the command-line interface (CLI) for a genome contig
classification pipeline. It orchestrates loading input data, performing
classification, generating summary statistics, and exporting results. 

The pipeline integrates contig size data with BLAST alignment results against
multiple reference databases to assign each contig to a biological category.

Workflow:
1. Parse command-line arguments
2. Load contig sizes and BLAST result files
3. Validate consistency of contig IDs across inputs
4. Classify contigs based on alignment coverage and size thresholds
5. Compute summary statistics per classification bin
6. Print a formatted summary report
7. Save the classification table (TSV)
8. Generate and save summary plots (PNG)

Parameters (via command-line arguments)
--------------------------------------
--sizes : str
    Path to contig sizes file (.txt). Must include contig identifiers and sizes in base pairs.

--mito : str
    Path to BLAST results against mitochondrial database (.tsv).

--apicomplexa : str
    Path to BLAST results against Apicomplexa database (.tsv).

--hepatozoon : str
    Path to BLAST results against Hepatozoon database (.tsv).
    These hits are merged into the Apicomplexa classification bin.

--sexual : str
    Path to BLAST results against sexual chromosome database (.tsv).

--outdir : str, optional
    Output directory for generated files (default: "output/").
    Created if it does not exist.

--min-coverage : float, optional
    Minimum alignment coverage fraction required to consider a BLAST hit
    valid (default: MIN_COVERAGE imported from contig_classifier).

--min-size : int, optional
    Minimum contig size (in base pairs) required for inclusion in the
    analysis (default: MIN_SIZE imported from contig_classifier).

Returns
-------
None
    This script is executed for its side effects. It does not return values
"""


import argparse
import pandas as pd

from contig_classifier import (
    load_sizes,
    load_blast,
    check_shared_ids,
    classify_contigs,
    MIN_COVERAGE,
    MIN_SIZE,
)
from summary import (
    compute_summary,
    save_charts,
    save_classification_table,
    print_summary,
)


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments for the contig classification pipeline.

    Parameters
    ----------
    None

    Ensures
    -------
    None

    Returns
    -------
    argparse.Namespace
        Parsed argument object with attributes ``sizes``, ``mito``,
        ``apicomplexa``, ``hepatozoon``, ``sexual``, ``outdir``,
        ``min_coverage``, and ``min_size``.
    """
    parser = argparse.ArgumentParser(
        description="Bin contigs from a draft genome assembly into biological categories."
    )
    parser.add_argument("--sizes",       required=True, help="Path to contig sizes file (.txt)")
    parser.add_argument("--mito",        required=True, help="BLAST results vs mitochondrial DB (.tsv)")
    parser.add_argument("--apicomplexa", required=True, help="BLAST results vs Apicomplexa DB (.tsv)")
    parser.add_argument("--hepatozoon",  required=True, help="BLAST results vs Hepatozoon DB (.tsv)")
    parser.add_argument("--sexual",      required=True, help="BLAST results vs sexual chromosomes DB (.tsv)")
    parser.add_argument("--outdir",      default="output", help="Output directory (default: output/)")
    parser.add_argument(
        "--min-coverage", type=float, default=MIN_COVERAGE,
        help=f"Minimum alignment coverage fraction (default: {MIN_COVERAGE})",
    )
    parser.add_argument(
        "--min-size", type=int, default=MIN_SIZE,
        help=f"Minimum contig size in bp to include (default: {MIN_SIZE})",
    )
    return parser.parse_args()


def load_all_blast(args: argparse.Namespace) -> dict[str, pd.DataFrame]:
    """
    Load all four BLAST result files and merge Hepatozoon hits into the Apicomplexa bin.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments. Must have attributes ``mito``,
        ``apicomplexa``, ``hepatozoon``, and ``sexual`` pointing to
        valid BLAST TSV files.

    Ensures
    -------
    None

    Returns
    -------
    dict[str, pd.DataFrame]
        Mapping of bin label to BLAST DataFrame with keys
        ``"Mitochondrion"``, ``"Apicomplexa"``, and ``"Sexual Chromosome"``.
        Hepatozoon hits are concatenated into the ``"Apicomplexa"`` entry.
    """
    mito_df = load_blast(args.mito,        "Mitochondrion")
    api_df  = load_blast(args.apicomplexa, "Apicomplexa")
    hep_df  = load_blast(args.hepatozoon,  "Apicomplexa")
    sex_df  = load_blast(args.sexual,      "Sexual Chromosome")
    blast_dfs = {
        "Mitochondrion":    mito_df,
        "Apicomplexa":      pd.concat([api_df, hep_df], ignore_index=True),
        "Sexual Chromosome": sex_df,
    }
    return blast_dfs


def run_pipeline(args: argparse.Namespace) -> None:
    """
    Execute the full contig classification pipeline end-to-end.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments containing all input file paths,
        output directory, and filtering thresholds.

    Ensures
    -------
    None

    Returns
    -------
    None
    """
    sizes_df  = load_sizes(args.sizes)
    blast_dfs = load_all_blast(args)
    check_shared_ids(blast_dfs, sizes_df)

    classified = classify_contigs(
        sizes_df,
        blast_dfs,
        min_coverage=args.min_coverage,
        min_size=args.min_size,
    )

    summary = compute_summary(classified)
    print_summary(summary)
    save_classification_table(classified, args.outdir)
    save_charts(summary, args.outdir)
    return


def main() -> None:
    """
    Entry point for the contigshw classification pipeline.

    Parameters
    ----------
    None

    Ensures
    -------
    None

    Returns
    -------
    None
    """
    args = parse_args()
    run_pipeline(args)
    return


if __name__ == "__main__":
    main()
