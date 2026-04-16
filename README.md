
# contigshw

A Python pipeline that classifies contigs from a draft genome assembly into biological bins — Mitochondrion, Apicomplexa, Sexual Chromosome, Diploid Chromosome, or Unclassified — using BLAST tabular output against multiple reference databases.

## Setup

1. Clone the repository:
```
git clone https://github.com/edecocke-uncc/contigshw.git
```
2. Go into the project folder:
```
cd contigshw
```
3. Create the environment:
```bash
conda env create -f environment.yml
```
4. Activate the environment:
```
conda activate contigshw
```

## Run

```
python3 src/main.py \
  --sizes data/Binsularis_contig_sizes.txt
  --mito data/Binsularis_BLAST_Mitochondrion.tsv
  --apicomplexa data/Binsularis_BLAST_Apicomplexa.tsv
  --hepatozoon data/Binsularis_BLAST_Hepatozoon.tsv
  --sexual data/Binsularis_BLAST_SexualChromosomes.tsv
  --outdir output/
```

All `--blast` inputs must be BLAST outfmt 6 TSV files with exactly 14 columns. The sizes file is a two-column whitespace-separated file with no header: contig name followed by length in base pairs.

### Optional flags

| Flag | Default | Description |
|---|---|---|
| `--min-coverage` | `0.01` | Minimum alignment coverage fraction (0–1) |
| `--min-size` | `3000` | Minimum contig length in bp |
| `--outdir` | `output/` | Output directory |

## What it does

* Loads a contig sizes file and validates it for duplicates, nulls, and non-positive lengths
* Loads four BLAST result files (mitochondrion, Apicomplexa, Hepatozoon, and sexual chromosome) and validates each against outfmt 6 requirements
* Merges Hepatozoon hits into the Apicomplexa bin, since *Hepatozoon* is a genus within phylum Apicomplexa; when a contig has hits in both tables, the hit with the highest bitscore drives classification
* Filters out contigs shorter than `--min-size` (default 3,000 bp) before classification
* Computes per-contig alignment coverage by merging overlapping HSPs (high-scoring segment pairs) before summing aligned bases, preventing inflation from overlapping alignments
* Applies the `--min-coverage` threshold to remove spurious short hits
* Resolves contigs with hits in multiple bins using a fixed priority order: Mitochondrion → Apicomplexa → Sexual Chromosome → Diploid Chromosome; within a bin, the highest-bitscore hit wins
* Labels any contig with no qualifying hit as `Unclassified`
* Writes a per-contig classification TSV and a two-panel summary bar chart

bins:
| Bin                     | Meaning                                      |
|-------------------------|----------------------------------------------|
| Mitochondrion           | Putative mitochondrial DNA                   |
| Apicomplexa             | Putative Apicomplexa parasite sequence       |
| Sexual Chromosome       | Putative chromosome W                        |
| Diploid Chromosome      | Putative autosomal (nuclear, non-W) sequence |

## Design decisions

### Coverage threshold (default: 1 %)

The minimum alignment coverage is set to **1 %** (`--min-coverage 0.01`). This balances two competing concerns:

**Sensitivity.** Organellar and parasite DNA in a draft assembly can be fragmented and interspersed with nuclear sequence. A stringent threshold (e.g., 50 %) would miss genuine mitochondrial or Apicomplexa contigs where only a portion of the contig carries the diagnostic sequence — for example, a nuclear-mitochondrial insertion or a chimeric contig produced by the assembler.

**Specificity.** Very short spurious hits (a few dozen bases on a multi-megabase contig) should not drive classification. At 1 %, a 100 kb contig needs at least 1 kb of aligned sequence; a 1 Mb contig needs at least 10 kb. This removes noise from incidental short matches while preserving biologically meaningful signal.

A different threshold (e.g., 0.1 %, 5 %, or a bitscore-based filter) may be appropriate depending on the assembly and reference databases; the value is easily overridden with `--min-coverage`.

### Hepatozoon integration

*Hepatozoon* is a genus within the phylum Apicomplexa. Hits from the Hepatozoon-specific database are therefore labeled with the same bin (`Apicomplexa`) as hits from the broader Apicomplexa database. When a contig has hits in both tables, the one with the highest bitscore is retained. This avoids double-counting and lets the more informative hit drive classification.

### Diploid Chromosome bin

Contigs classified as "Diploid Chromosome" would be those that appear only in the Sexual Chromosome BLAST table but fail a criterion distinguishing autosomes from the W chromosome. However, the current BLAST data only includes a chromosome W database — there is no autosome-specific database to provide positive evidence for the diploid bin. As a result, no contigs are assigned to "Diploid Chromosome" in this run. Contigs without hits in any database are labeled "Unclassified"; many of these are likely autosomal.

A more refined approach would involve additional databases (e.g., known autosomal sequences from a closely related species) or k-mer-based coverage analysis to separate W-linked from autosomal contigs.

### Overlap merging

When computing alignment coverage, overlapping HSPs on the query are merged into non-overlapping intervals before summing. This prevents inflating coverage when BLAST returns multiple partially overlapping alignments to different subjects in the same database.

## Output

Two files are written to `--outdir` (default: `output/`) after each run:

* `classification.tsv` — one row per eligible contig with columns `qseqid`, `size_bp`, `bin`, `sseqid`, `bitscore`, `evalue`, `pident`, and `coverage`
* `summary_charts.png` — side-by-side bar charts showing contig count and total assembly size (Mbp) per bin
