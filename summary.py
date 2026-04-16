import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

BINS = ["Mitochondrion", "Apicomplexa", "Sexual Chromosome", "Diploid Chromosome", "Unclassified"]
COLORS = ["#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#999999"]


def compute_summary(df):
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


def annotate_bars(ax, bars, values, fmt):
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


def plot_contig_counts(summary, ax):
    bars = ax.bar(
        summary["Bin"], summary["n_contigs"],
        color=COLORS, edgecolor="white", linewidth=0.8,
    )
    ax.set_title("Contigs per bin", fontsize=13, fontweight="bold")
    ax.set_ylabel("Number of contigs")
    ax.set_xticklabels(summary["Bin"], rotation=30, ha="right")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    annotate_bars(ax, bars, summary["n_contigs"], "{:,}")
    return


def plot_total_bp(summary, ax):
    bars = ax.bar(
        summary["Bin"], summary["total_mbp"],
        color=COLORS, edgecolor="white", linewidth=0.8,
    )
    ax.set_title("Total assembly size per bin", fontsize=13, fontweight="bold")
    ax.set_ylabel("Total base pairs (Mbp)")
    ax.set_xticklabels(summary["Bin"], rotation=30, ha="right")
    annotate_bars(ax, bars, summary["total_mbp"], "{:.1f}")
    return


def save_charts(summary, output_dir):
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


def save_classification_table(df, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "classification.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"[INFO] Classification table saved to {out_path}")
    return out_path


def print_summary(summary):
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
