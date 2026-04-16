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


def parse_args():
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


def load_all_blast(args):
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


def run_pipeline(args):
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


def main():
    args = parse_args()
    run_pipeline(args)
    return


if __name__ == "__main__":
    main()
