import os
import sys
import pandas as pd
import numpy as np
from typing import Optional


BLAST_COLS = [
    "qseqid", "staxids", "bitscore", "qseqid2", "sseqid",
    "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue",
]

PRIORITY_ORDER = ["Mitochondrion", "Apicomplexa", "Sexual Chromosome", "Diploid Chromosome"]
MIN_SIZE = 3000
MIN_COVERAGE = 0.01


def validate_file_exists(path: str) -> None:
    """
    Validate that a file path exists, is readable, and is non-empty.

    Parameters
    ----------
    path : str
        Absolute or relative path to the file to validate.
        Prints a warning if a BLAST file does not have a .tsv extension.

    Ensures
    -------
    - The file at ``path`` exists on disk.
    - The file at ``path`` is readable by the current process.
    - The file at ``path`` has a non-zero size.

    Returns
    -------
    None
    """
    if not path.endswith(".tsv") and "BLAST" in os.path.basename(path):
        print(f"[WARNING] Expected .tsv extension for BLAST file: {path}")
    if not os.path.exists(path):
        sys.exit(f"[ERROR] File does not exist: {path}")
    if not os.access(path, os.R_OK):
        sys.exit(f"[ERROR] File is not readable: {path}")
    if os.path.getsize(path) == 0:
        sys.exit(f"[ERROR] File is empty: {path}")
    return


def validate_outfmt6(df: pd.DataFrame, path: str) -> None:
    """
    Validate that a DataFrame conforms to BLAST tabular outfmt 6 expectations.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame loaded from a BLAST outfmt 6 TSV file.
    path : str
        File path used in error messages to identify the source file.

    Ensures
    -------
    - ``df`` has exactly 14 columns.
    - All required outfmt 6 column names are present.
    - The ``qseqid`` column contains no null values.
    - The ``sseqid`` column contains no null values.
    - The ``bitscore`` column is not entirely null.
    - The ``evalue`` column is not entirely null.
    - The ``length`` column is not entirely null.

    Returns
    -------
    None
    """
    if df.shape[1] != 14:
        sys.exit(
            f"[ERROR] {path} has {df.shape[1]} columns; outfmt6 requires exactly 14."
        )
    required = ["qseqid", "sseqid", "bitscore", "evalue", "length", "pident", "qstart", "qend"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"[ERROR] {path} is missing required outfmt6 columns: {missing}")
    if df["qseqid"].isnull().any():
        sys.exit(f"[ERROR] Null query identifiers (qseqid) found in {path}")
    if df["sseqid"].isnull().any():
        sys.exit(f"[ERROR] Null subject identifiers (sseqid) found in {path}")
    if df["bitscore"].isnull().all():
        sys.exit(f"[ERROR] All bitscore values are null in {path}")
    if df["evalue"].isnull().all():
        sys.exit(f"[ERROR] All evalue values are null in {path}")
    if df["length"].isnull().all():
        sys.exit(f"[ERROR] All alignment length values are null in {path}")
    return


def check_shared_ids(blast_dfs: dict[str, pd.DataFrame], sizes_df: pd.DataFrame) -> None:
    """
    Verify that each BLAST result DataFrame shares at least one contig ID with the sizes file.

    Parameters
    ----------
    blast_dfs : dict[str, pd.DataFrame]
        Mapping of bin label to BLAST result DataFrame. Each DataFrame must
        contain a ``qseqid`` column of contig identifiers.
    sizes_df : pd.DataFrame
        Contig sizes DataFrame containing a ``contig_name`` column.

    Ensures
    -------
    - Every BLAST DataFrame in ``blast_dfs`` has at least one contig ID that
      also appears in ``sizes_df["contig_name"]``.

    Returns
    -------
    None
    """
    size_ids = set(sizes_df["contig_name"])
    for label, df in blast_dfs.items():
        blast_ids = set(df["qseqid"])
        shared = blast_ids & size_ids
        if len(shared) == 0:
            sys.exit(
                f"[ERROR] No shared contig IDs between '{label}' BLAST results and sizes file. "
                "Check that the correct files were provided."
            )
        not_found = blast_ids - size_ids
        print(f"[QC] {label}: {len(shared)} shared IDs | {len(not_found)} query IDs not in sizes file.")
    return


def validate_sizes(df: pd.DataFrame) -> None:
    """
    Validate the contig sizes DataFrame for integrity issues.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns ``contig_name`` (str) and ``size_bp`` (int).

    Ensures
    -------
    - No null values exist in the ``contig_name`` column.
    - No duplicate contig names exist in the ``contig_name`` column.
    - No null values exist in the ``size_bp`` column.
    - All values in ``size_bp`` are strictly positive.

    Returns
    -------
    None
    """
    if df["contig_name"].isnull().any():
        sys.exit("[ERROR] Null contig identifiers found in sizes file.")
    if df["contig_name"].duplicated().any():
        dupes = df[df["contig_name"].duplicated()]["contig_name"].tolist()
        sys.exit(f"[ERROR] Duplicate contig names in sizes file: {dupes[:5]}")
    if df["size_bp"].isnull().any():
        sys.exit("[ERROR] Null size values found in sizes file.")
    if (df["size_bp"] <= 0).any():
        sys.exit("[ERROR] Non-positive contig sizes found in sizes file.")
    return


def load_sizes(path: str) -> pd.DataFrame:
    """
    Load and validate a whitespace-delimited contig sizes file.

    Parameters
    ----------
    path : str
        Path to a two-column text file with contig name and size in base pairs.
        Columns are expected to be whitespace-separated with no header.

    Ensures
    -------
    - The file passes ``validate_file_exists`` checks.
    - The loaded DataFrame passes ``validate_sizes`` checks.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ``contig_name`` (str) and ``size_bp`` (int).
    """
    validate_file_exists(path)
    df = pd.read_csv(path, sep=r"\s+", header=None, names=["contig_name", "size_bp"])
    validate_sizes(df)
    print(f"[QC] Sizes file: {len(df)} contigs loaded from {path}")
    return df


def load_blast(path: str, label: str) -> pd.DataFrame:
    """
    Load and validate a BLAST outfmt 6 TSV file and attach a bin label.

    Parameters
    ----------
    path : str
        Path to the BLAST tabular output file (outfmt 6, 14 columns).
    label : str
        Bin label to assign to all hits in this file (e.g. ``"Mitochondrion"``).

    Ensures
    -------
    - The file passes ``validate_file_exists`` checks.
    - The loaded DataFrame passes ``validate_outfmt6`` checks.
    - ``qstart`` and ``qend`` are reoriented so that ``qstart <= qend``
      for every row (handles reverse-strand hits).

    Returns
    -------
    pd.DataFrame
        Cleaned BLAST DataFrame with a ``source`` column set to ``label``
        and the redundant ``qseqid2`` column removed.
    """
    validate_file_exists(path)
    df = pd.read_csv(path, sep="\t", header=None, names=BLAST_COLS)
    validate_outfmt6(df, path)
    df = df.drop(columns=["qseqid2"])
    df["qstart"], df["qend"] = (
        np.minimum(df["qstart"], df["qend"]),
        np.maximum(df["qstart"], df["qend"]),
    )
    df["source"] = label
    print(f"[QC] {label}: {len(df)} BLAST hits loaded from {path}")
    return df


def merge_intervals(intervals: list[tuple[int, int]]) -> list[list[int]]:
    """
    Merge a list of possibly overlapping integer intervals into non-overlapping ones.

    Parameters
    ----------
    intervals : list[tuple[int, int]]
        List of ``(start, end)`` coordinate pairs (both inclusive, 1-based).
        May be unsorted or contain overlapping entries.

    Ensures
    -------
    None

    Returns
    -------
    list[list[int]]
        Sorted list of non-overlapping ``[start, end]`` intervals that cover
        exactly the same positions as the input. Returns an empty list when
        ``intervals`` is empty.
    """
    if not intervals:
        return []
    sorted_iv = sorted(intervals, key=lambda x: x[0])
    merged = [list(sorted_iv[0])]
    for start, end in sorted_iv[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def coverage_for_group(group: pd.DataFrame) -> float:
    """
    Compute the fraction of a contig covered by merged BLAST HSPs.

    Parameters
    ----------
    group : pd.DataFrame
        Subset of a merged BLAST + sizes DataFrame for a single contig.
        Must contain columns ``size_bp``, ``qstart``, and ``qend``.

    Ensures
    -------
    None

    Returns
    -------
    float
        Coverage fraction in the range ``[0.0, 1.0]``. Returns ``0.0`` if
        the contig size is zero to avoid division by zero.
    """
    size = group["size_bp"].iloc[0]
    intervals = list(zip(group["qstart"].tolist(), group["qend"].tolist()))
    merged = merge_intervals(intervals)
    aligned_bp = sum(e - s + 1 for s, e in merged)
    return float(aligned_bp) / size if size > 0 else 0.0


def compute_all_coverages(blast_df: pd.DataFrame, sizes_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-contig alignment coverage for all contigs in a BLAST result.

    Parameters
    ----------
    blast_df : pd.DataFrame
        BLAST result DataFrame containing columns ``qseqid``, ``qstart``,
        and ``qend``.
    sizes_df : pd.DataFrame
        Contig sizes DataFrame with columns ``contig_name`` and ``size_bp``.

    Ensures
    -------
    None

    Returns
    -------
    pd.DataFrame
        Two-column DataFrame with ``qseqid`` and ``coverage`` (float fraction).
        Only contigs present in both ``blast_df`` and ``sizes_df`` are included.
    """
    merged = blast_df.merge(sizes_df, left_on="qseqid", right_on="contig_name", how="inner")
    cov_series = merged.groupby("qseqid", group_keys=False).apply(lambda g: coverage_for_group(g), include_groups=False)
    cov_df = cov_series.reset_index()
    cov_df.columns = ["qseqid", "coverage"]
    return cov_df


def get_best_hit_per_contig(df: pd.DataFrame) -> pd.DataFrame:
    """
    Retain only the highest-bitscore BLAST hit for each contig.

    Parameters
    ----------
    df : pd.DataFrame
        BLAST result DataFrame containing at least ``qseqid`` and ``bitscore``
        columns.

    Ensures
    -------
    None

    Returns
    -------
    pd.DataFrame
        Deduplicated DataFrame with one row per unique ``qseqid``, keeping the
        row with the highest ``bitscore``. Row index is reset.
    """
    return (
        df.sort_values("bitscore", ascending=False)
        .drop_duplicates(subset="qseqid")
        .reset_index(drop=True)
    )


def build_candidate_table(
    blast_dfs: dict[str, pd.DataFrame],
    sizes_df: pd.DataFrame,
    min_coverage: float,
) -> pd.DataFrame:
    """
    Build a combined table of all BLAST hits that meet the coverage threshold.

    Parameters
    ----------
    blast_dfs : dict[str, pd.DataFrame]
        Mapping of bin label to BLAST result DataFrame.
    sizes_df : pd.DataFrame
        Contig sizes DataFrame with columns ``contig_name`` and ``size_bp``.
    min_coverage : float
        Minimum alignment coverage fraction (between 0 and 1) required for a
        contig to be included as a candidate.

    Ensures
    -------
    None

    Returns
    -------
    pd.DataFrame
        Combined DataFrame of all hits that pass the coverage filter, with
        columns ``qseqid``, ``bin``, ``sseqid``, ``bitscore``, ``evalue``,
        ``pident``, and ``coverage``. Returns an empty DataFrame with those
        columns if no hits pass the threshold.
    """
    frames = []
    for label, df in blast_dfs.items():
        cov_df = compute_all_coverages(df, sizes_df)
        df2 = df.merge(cov_df, on="qseqid")
        df2 = df2[df2["coverage"] >= min_coverage].copy()
        df2["bin"] = label
        frames.append(df2)
    if not frames:
        empty_cols = ["qseqid", "bin", "sseqid", "bitscore", "evalue", "pident", "coverage"]
        return pd.DataFrame(columns=empty_cols)
    return pd.concat(frames, ignore_index=True)


def apply_priority_and_best_hit(candidates: pd.DataFrame) -> pd.DataFrame:
    """
    Resolve multi-bin conflicts by applying the bin priority order and best bitscore.

    Parameters
    ----------
    candidates : pd.DataFrame
        Combined candidate table from ``build_candidate_table``. Must contain
        columns ``qseqid``, ``bin``, and ``bitscore``.

    Ensures
    -------
    None

    Returns
    -------
    pd.DataFrame
        One row per unique ``qseqid`` with columns ``qseqid``, ``bin``,
        ``sseqid``, ``bitscore``, ``evalue``, ``pident``, and ``coverage``.
        When a contig has hits in multiple bins, the bin earliest in
        ``PRIORITY_ORDER`` is preferred; ties within a bin are broken by
        highest bitscore.
    """
    priority_map = {b: i for i, b in enumerate(PRIORITY_ORDER)}
    candidates = candidates.copy()
    candidates["priority"] = candidates["bin"].map(priority_map)
    candidates = candidates.sort_values(
        ["qseqid", "priority", "bitscore"],
        ascending=[True, True, False],
    )
    best = candidates.drop_duplicates(subset="qseqid", keep="first")
    keep_cols = ["qseqid", "bin", "sseqid", "bitscore", "evalue", "pident", "coverage"]
    return best[keep_cols].reset_index(drop=True)


def classify_contigs(
    sizes_df: pd.DataFrame,
    blast_dfs: dict[str, pd.DataFrame],
    min_coverage: float = MIN_COVERAGE,
    min_size: int = MIN_SIZE,
) -> pd.DataFrame:
    """
    Classify all contigs into biological bins using BLAST evidence.

    Parameters
    ----------
    sizes_df : pd.DataFrame
        Contig sizes DataFrame with columns ``contig_name`` and ``size_bp``.
    blast_dfs : dict[str, pd.DataFrame]
        Mapping of bin label to BLAST result DataFrame.
    min_coverage : float, optional
        Minimum alignment coverage fraction required to consider a hit
        (default: ``MIN_COVERAGE`` = 0.01).
    min_size : int, optional
        Minimum contig length in base pairs to include in classification
        (default: ``MIN_SIZE`` = 3000).

    Ensures
    -------
    - ``min_coverage`` is between 0 and 1 (inclusive).
    - ``min_size`` is a positive integer.

    Returns
    -------
    pd.DataFrame
        One row per contig that meets the size threshold, with columns
        ``qseqid``, ``size_bp``, ``bin``, ``sseqid``, ``bitscore``,
        ``evalue``, ``pident``, and ``coverage``. Contigs with no qualifying
        BLAST hits are labeled ``"Unclassified"`` with NaN for hit columns.
    """
    if min_coverage < 0 or min_coverage > 1:
        sys.exit(f"[ERROR] --min-coverage must be between 0 and 1, got {min_coverage}")
    if min_size <= 0:
        sys.exit(f"[ERROR] --min-size must be a positive integer, got {min_size}")

    eligible = sizes_df[sizes_df["size_bp"] >= min_size].copy()
    print(f"[INFO] {len(eligible)} contigs >= {min_size} bp retained for classification.")

    candidates = build_candidate_table(blast_dfs, eligible, min_coverage)
    base = eligible[["contig_name", "size_bp"]].rename(columns={"contig_name": "qseqid"})

    if candidates.empty:
        base["bin"] = "Unclassified"
        for col in ["sseqid", "bitscore", "evalue", "pident", "coverage"]:
            base[col] = np.nan
        return base

    best = apply_priority_and_best_hit(candidates)
    final = base.merge(best, on="qseqid", how="left")
    final["bin"] = final["bin"].fillna("Unclassified")
    return final
