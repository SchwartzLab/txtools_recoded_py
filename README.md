# txtools (Python)

A fast, easy-to-install **Python reimplementation** of the
[txtools](https://github.com/AngelCampos/txtools) R package
(Garcia-Campos *et al.*, *Nucleic Acids Research* 2024).

txtools converts genomic RNA-seq alignments into **transcriptomic coordinates**
and summarizes them into a per-nucleotide table capturing **coverage**,
**read-starts / read-ends**, **nucleotide frequencies**, **deletions**, and
**inserts** — the readouts used for epitranscriptomic (RNA modification),
RNA-structure, and RNA:protein-interaction analyses.

This rewrite keeps the **function names and workflow you already know**, but
replaces the R/Bioconductor stack with `pysam` + `numpy` + `pandas` for speed
and a clean install (no version clashes), and ships a single `txtools`
command-line tool.

> Status: core pipeline + metrics. Plotting and statistical tests
> (`tx_plot_*`, `tx_test_*`) are intentionally out of scope for this version —
> the TSV table feeds directly into your own plotting/stats in Python or R.

---

## Install (conda-first — recommended)

Conda isolates dependencies per-tool, which is exactly what avoids the
version-clash pain of the R package. `pysam` installs precompiled from
bioconda (no compiler needed).

```bash
conda env create -f conda/environment.yml   # creates env "txtools"
conda activate txtools
# the env file already does an editable install; otherwise:
pip install -e .
```

Or with an existing environment that has `pysam`, `numpy`, `pandas`:

```bash
pip install -e .          # exposes the `txtools` command + Python API
```

### Standalone binary (optional, on-demand)

For machines without conda you can build a single self-contained executable
(platform-specific, ~hundreds of MB; rebuild per fix):

```bash
pip install -e ".[binary]"
make binary               # -> dist/txtools
```

---

## Command-line usage

One shot: BAM + BED12 + FASTA → transcriptomic TSV.

```bash
txtools makeDT \
    --bam   reads.bam \            # coordinate-sorted, indexed (.bai auto-created)
    --bed   genes.bed12 \          # gene models (BED-12 keeps exon structure; BED-6 ok)
    --fasta genome.fa \            # reference (adds the refSeq column)
    --table covNucFreq \           # coverage | nucFreq | covNucFreq
    --paired \                     # or --single
    --strand-mode 1 \              # 1 = directional (default); 2 = dUTP/TruSeq
    --min-reads 50 \
    --threads 8 \
    --add misincRate \             # derived metrics, repeat --add per metric
    --add startRatio \
    --add geneRegion \
    --add misincRateNucSpec:C,T \
    --add motifPresence:DRACH,center \
    -o out.tsv                     # ".tsv.gz" is compressed; "-" = stdout
```

### `--add` metrics

Repeat `--add` once per metric. Metrics taking arguments use `name:arg1,arg2`:

| Metric | Form |
| --- | --- |
| `startRatio`, `endRatio`, `*1bpDS`, `*1bpUS` | `--add startRatio` |
| `nucTotal`, `misincCount`, `misincRate` | `--add misincRate` |
| `geneRegion`, `relTxPos`, `pos`, `refSeqDT` | `--add geneRegion` |
| `misincRateNucSpec` | `--add misincRateNucSpec:C,T` |
| `motifPresence` | `--add motifPresence:DRACH,center` (or `:DRACH,3` / `:DRACH,all`) |
| `rollingMean` | `--add rollingMean:cov,50,center` |

---

## Python API

The functions mirror the R package, so existing mental models carry over:

```python
import txtools

ga     = txtools.tx_load_bed("genes.bed12")
genome = txtools.tx_load_genome("genome.fa")
bam    = txtools.tx_load_bam("reads.bam", paired_end=True, load_seq=True)

reads  = txtools.tx_reads(bam, ga, min_reads=50, with_seq=True, n_cores=8)
DT     = txtools.tx_makeDT_covNucFreq(reads, ga, genome=genome)

DT = (DT
      .pipe(txtools.tx_add_misincRate, min_nuc_reads=20)
      .pipe(txtools.tx_add_startRatio)
      .pipe(txtools.tx_add_geneRegion, ga)
      .pipe(txtools.tx_add_motifPresence, "DRACH", nuc_positions="center"))

DT.to_csv("out.tsv", sep="\t", index=False, na_rep="NA")
```

The txDT is a `pandas.DataFrame` (one row per transcriptomic position):

| column | meaning |
| --- | --- |
| `chr`, `gencoor`, `strand`, `gene`, `txcoor` | genomic + transcriptomic coordinates (1-based) |
| `refSeq` | reference base (strand-aware), if a genome was given |
| `cov`, `start_5p`, `end_3p` | insert coverage, read-start and read-end counts |
| `A C G T [M R W S Y K] N - .` | nucleotide-frequency pileup (`-` deletion, `.` insert) |

---

## Differences from the R package (by design)

* **Inputs** stay BAM (coordinate-sorted + indexed), BED-12/6, FASTA. `tx_reads`
  fetches reads per gene via the BAM index instead of loading the whole file
  into memory.
* **`strandMode = 2`**: transcript boundaries are derived from the fragment's
  genomic extent + gene strand. This is **numerically identical to R for
  `strandMode = 1`** (the common directional case) and **correct for mode 2**,
  where the original mate-role formula could invert coordinates.
* **Splice-junction check** is applied per-read (only a read whose junctions
  disagree with the annotation is dropped), rather than the R code's
  all-or-nothing drop of every spliced read in a gene.
* **Overlap consensus** of mate pairs reproduces the well-defined IUPAC cases;
  rare discordant overlaps involving `N` resolve to `N`.
* Out of scope (for now): `tx_plot_*`, `tx_test_*`, and the unassigned-read
  dump.

---

## Development

```bash
pip install -e ".[dev]"
make test          # pytest on a self-contained toy dataset
```

The test suite builds a tiny genome/annotation/BAM on the fly and checks
coordinate conversion (+ and − strand), spliced reads, insert coverage,
mismatch pileups, IUPAC splitting, and every metric.
```
