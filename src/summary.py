#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Erin Nicole Decocker
# edecocke@charlotte.edu
# ID: 801442694
# AI usage acknowledgment in READme file.

"""
summary.py

This script provides utilities for summarizing, visualizing, and exporting
genome contig classification results. It operates on a per-contig
classification DataFrame (e.g., output from a `classify_contigs` function)
and produces both tabular summaries and publication-ready plots.

Core functionality includes:
- Aggregating contigs into predefined biological bins
- Computing counts and total assembly size per bin
- Generating bar charts for contig counts and total base pairs (in Mbp)
- Saving summary plots as a PNG file
- Exporting the full classification table as a TSV file
- Printing a formatted summary report to standard output

Expected input:
- A pandas DataFrame containing at minimum:
    - `qseqid` : contig identifier
    - `bin`    : assigned classification bin
    - `size_bp`: contig size in base pairs

Output:
- Summary statistics per bin (counts and total Mbp)
- Visualization file: `summary_charts.png`
- Classification table: `classification.tsv`
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.axes import Axes
from matplotlib.container import BarContainer


BINS = ["Mitochondrion", "Apicomplexa", "Sexual Chromosome", "Diploid Chromosome", "Unclassified"]
COLORS = ["#D95F02", "#1B9E77", "#00CECE", "#E7298A", "#CE0000"]


def compute_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate classification results into per-bin contig counts and total base pairs.

    Parameters
    ----------
    df : pd.DataFrame
        Classification DataFrame output by ``classify_contigs``. Must contain
        columns ``bin`` and ``size_bp`` with one row per contig, and ``qseqid``
        as the contig identifier.

    Ensures
    -------
    None

    Returns
    -------
    pd.DataFrame
        Summary DataFrame with one row per bin (ordered by ``BINS``), containing
        columns ``Bin`` (str), ``n_contigs`` (int), ``total_bp`` (int), and
        ``total_mbp`` (float). Bins with no contigs are included with zero counts.
    """
    summary = (
        df.groupby("bin", observed=False)
        .agg(
            n_contigs=("qseqid", "count"),
            total_bp=("size_bp", "sum"),
        )
        .reindex(BINS, fill_value=0)
        .reset_index()
    )
    summary.rename(columns={"bin": "Bin"}, inplace=True)
    summary["total_mbp"] = summary["total_bp"] / 1e6
    return summary


def annotate_bars(ax: Axes, bars: BarContainer, values: pd.Series, fmt: str) -> None:
    """
    Place a text label above each bar in a bar chart.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The Axes object containing the bar chart.
    bars : matplotlib.container.BarContainer
        The collection of bar patches returned by ``ax.bar()``.
    values : pd.Series
        Numeric values to format and display, one per bar.
    fmt : str
        Python format string applied to each value (e.g. ``"{:,}"`` or
        ``"{:.1f}"``).

    Ensures
    -------
    None

    Returns
    -------
    None
    """
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + ax.get_ylim()[1] * 0.01,
            fmt.format(val),
            ha="center",
            va="bottom",
            fontsize=9,
        )
    return


def plot_contig_counts(summary: pd.DataFrame, ax: Axes) -> None:
    """
    Draw a bar chart of contig counts per bin on a given Axes.

    Parameters
    ----------
    summary : pd.DataFrame
        Summary DataFrame from ``compute_summary``. Must contain columns
        ``Bin`` and ``n_contigs``.
    ax : matplotlib.axes.Axes
        The Axes object on which to draw the chart.

    Ensures
    -------
    None

    Returns
    -------
    None
    """
    bars = ax.bar(
        summary["Bin"], summary["n_contigs"],
        color=COLORS, edgecolor="white", linewidth=0.8,
    )
    ax.set_title("Contigs per Bin", fontsize=13, fontweight="bold")
    ax.set_ylabel("Number of contigs")
    ax.set_xticks(range(len(summary["Bin"]))) 
    ax.set_xticklabels(summary["Bin"], rotation=30, ha="right")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    annotate_bars(ax, bars, summary["n_contigs"], "{:,}")
    return


def plot_total_bp(summary: pd.DataFrame, ax: Axes) -> None:
    """
    Draw a bar chart of total assembly size in megabase pairs per bin on a given Axes.

    Parameters
    ----------
    summary : pd.DataFrame
        Summary DataFrame from ``compute_summary``. Must contain columns
        ``Bin`` and ``total_mbp``.
    ax : matplotlib.axes.Axes
        The Axes object on which to draw the chart.

    Ensures
    -------
    None

    Returns
    -------
    None
    """
    bars = ax.bar(
        summary["Bin"], summary["total_mbp"],
        color=COLORS, edgecolor="white", linewidth=0.8,
    )
    ax.set_title("Assembly Size per Bin", fontsize=13, fontweight="bold")
    ax.set_ylabel("Base pairs (Mbp)")
    ax.set_xticks(range(len(summary["Bin"]))) 
    ax.set_xticklabels(summary["Bin"], rotation=30, ha="right")
    annotate_bars(ax, bars, summary["total_mbp"], "{:.1f}")
    return


def save_charts(summary: pd.DataFrame, output_dir: str) -> str:
    """
    Render and save a two-panel summary bar chart to disk.

    Parameters
    ----------
    summary : pd.DataFrame
        Summary DataFrame from ``compute_summary``.
    output_dir : str
        Directory path where the chart PNG will be written. Created if it
        does not already exist.

    Ensures
    -------
    None

    Returns
    -------
    str
        Absolute or relative path to the saved PNG file
        (``<output_dir>/summary_charts.png``).
    """
    os.makedirs(output_dir, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    plot_contig_counts(summary, ax1)
    plot_total_bp(summary, ax2)
    fig.tight_layout(pad=2.0)
    out_path = os.path.join(output_dir, "summary_charts.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Charts saved to {out_path}")
    return out_path


def save_classification_table(df: pd.DataFrame, output_dir: str) -> str:
    """
    Write the full per-contig classification table to a TSV file.

    Parameters
    ----------
    df : pd.DataFrame
        Classification DataFrame output by ``classify_contigs``.
    output_dir : str
        Directory path where the TSV will be written. Created if it does
        not already exist.

    Ensures
    -------
    None

    Returns
    -------
    str
        Absolute or relative path to the saved TSV file
        (``<output_dir>/classification.tsv``).
    """
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "classification.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"[INFO] Classification table saved to {out_path}")
    return out_path


def print_summary(summary: pd.DataFrame) -> None:
    """
    Print a formatted classification summary table to stdout.

    Parameters
    ----------
    summary : pd.DataFrame
        Summary DataFrame from ``compute_summary``. Must contain columns
        ``Bin``, ``n_contigs``, and ``total_mbp``.

    Ensures
    -------
    None

    Returns
    -------
    None
    """
    total = summary["n_contigs"].sum()
    classified = summary.loc[summary["Bin"] != "Unclassified", "n_contigs"].sum()
    unclassified_row = summary.loc[summary["Bin"] == "Unclassified", "n_contigs"]
    unclassified = int(unclassified_row.values[0]) if not unclassified_row.empty else 0

    print("\n========== CLASSIFICATION SUMMARY ==========")
    for _, row in summary.iterrows():
        print(f"  {row['Bin']:<25} {int(row['n_contigs']):>6} contigs   {row['total_mbp']:>10.1f} Mbp")
    print(f"\n  Total eligible contigs : {total}")
    print(f"  Classified             : {classified}")
    print(f"  Unclassified           : {unclassified}")
    print("============================================\n")
    return
