# Running txtools (step by step)

This guide is written for the Schwartz lab on WEXAC, but works anywhere.

---

## 0. What you need

Three inputs (all standard):

| Input | What | Notes |
|-------|------|-------|
| **BAM** | aligned reads | must be **coordinate-sorted** and **indexed** (`.bai`). A `.bai` is created automatically if missing. |
| **BED** | gene models | BED-12 (keeps exon structure) or BED-6 (single-exon). Gene names (column 4) must be unique. |
| **FASTA** | reference genome | only needed if you want the `refSeq` column and misincorporation metrics. A `.fai` is created automatically. |

Chromosome names must match across the three files (e.g. all `chr19`, not a mix of `19`/`chr19`).

---

## 1. Install (once)

### Option A — it's already installed in the `env` conda env
```bash
/home/labs/schwartzlab/schwartz/.conda/envs/env/bin/txtools --version
```
If that prints `txtools 0.1.0`, you can use that `txtools` directly (prefix with the path, or `conda activate env`).

### Option B — a clean, dedicated conda env (recommended, avoids clashes)
```bash
cd /home/labs/schwartzlab/schwartz/tools/txtools_recoded
conda env create -f conda/environment.yml      # creates env "txtools"
conda activate txtools
txtools --version
```

### Option C — pip into any env that has pysam/numpy/pandas
```bash
pip install -e /home/labs/schwartzlab/schwartz/tools/txtools_recoded
```

---

## 2. Run it (command line)

The one-shot command is `txtools makeDT`. Minimal example:

```bash
txtools makeDT \
    --bam   /path/to/reads.bam \
    --bed   /path/to/genes.bed12 \
    --fasta /path/to/genome.fa \
    -o      out.tsv
```

Full example (single-end STAMP, like the validation run):

```bash
txtools makeDT \
    --bam   sample_AlignedReads.bam \
    --bed   refGene.bed12 \
    --fasta genome.fa \
    --table covNucFreq \          # coverage | nucFreq | covNucFreq
    --single \                    # use --paired for paired-end (default)
    --strand-mode 2 \             # 1 = directional (default); 2 = dUTP/TruSeq/STAMP
    --min-reads 50 \
    --threads 8 \
    --add misincRate \            # repeat --add per metric
    --add misincRateNucSpec:C,T \ # STAMP C→T editing rate
    -o sample_txDT.tsv
```

Output is a TSV (one row per transcriptomic position). Use `-o out.tsv.gz` to gzip, or `-o -` for stdout.

### Picking `--strand-mode`
If unsure, run once each way with `--table coverage --min-reads 50` and keep whichever assigns more
reads (printed as "Output: N fragments in M gene models"). Mode **2** is right for STAMP / TruSeq-stranded.

### `--add` metric cheatsheet
```
--add startRatio            --add endRatio
--add misincRate            --add nucTotal
--add geneRegion            --add relTxPos
--add misincRateNucSpec:C,T          # ref,mis nucleotides
--add motifPresence:DRACH,center     # motif,position (center|all|<int>)
--add rollingMean:cov,50,center      # column,window,align
```

---

## 3. Run it (Python / notebooks)

```python
import txtools

ga     = txtools.tx_load_bed("genes.bed12")
genome = txtools.tx_load_genome("genome.fa")
bam    = txtools.tx_load_bam("reads.bam", paired_end=False, strand_mode=2, load_seq=True)

reads  = txtools.tx_reads(bam, ga, min_reads=50, with_seq=True, n_cores=8)
DT     = txtools.tx_makeDT_covNucFreq(reads, ga, genome=genome)        # a pandas DataFrame
DT     = txtools.tx_add_misincRate(DT)
DT     = txtools.tx_add_misincRateNucSpec(DT, ref_nuc="C", mis_nuc="T")

DT.to_csv("out.tsv", sep="\t", index=False, na_rep="NA")
```

---

## 4. Reproduce the chr19 validation (Python vs R)

Everything is in `validation/`:

```bash
V=/home/labs/schwartzlab/schwartz/tools/txtools_recoded/validation
BIN=/home/labs/schwartzlab/schwartz/.conda/envs/env/bin

# (a) Python tool
$BIN/txtools makeDT --bam $V/chr19.bam --bed $V/chr19_refGene.bed12 --fasta $V/chr19.fa \
    --table covNucFreq --single --strand-mode 2 --min-reads 50 --threads 8 \
    --add misincRate --add misincRateNucSpec:C,T -o $V/py_chr19.tsv

# (b) Reference R txtools 1.0.6  (R 4.3.1 + Bioconductor 3.17; txtools isolated in its own R_LIBS
#     so compiled deps come from the consistent site library)
export R_LIBS=$V/txlib
Rscript $V/run_R_txtools.R $V 8

# (c) Comparison HTML
$BIN/python $V/compare_report.py $V
#   -> $V/txtools_comparison_chr19.html
```

To regenerate the chr19 inputs from scratch (subset BAM, chr19 FASTA, refGene BED12) see the
commands recorded at the top of this project's history; the originals are already in `validation/`.

---

## Troubleshooting: `ModuleNotFoundError: No module named 'numpy.core._multiarray_umath'`

On WEXAC, the shell auto-loads easybuild modules (SciPy-bundle, GDAL, …) that inject a
**Python 3.10** `PYTHONPATH` and put easybuild's Python ahead of conda on `PATH`. Inside a conda
env on a *different* Python (e.g. 3.13) this makes Python import the wrong, incompatible numpy.

Fix (one line, before running):
```bash
unset PYTHONPATH
```
Made permanent for the `txtools` env via an activate hook
(`$CONDA_PREFIX/etc/conda/activate.d/zzz_unset_pythonpath.sh`) that runs `unset PYTHONPATH` and puts
the env first on `PATH`. So simply:
```bash
conda activate txtools     # hook auto-clears PYTHONPATH and fixes PATH
python -c "import numpy"    # -> .../envs/txtools/.../numpy   (the env's own)
```
If you create the env yourself elsewhere, add that one-line hook (or just `unset PYTHONPATH`).

## Tips
- Start on a subset (one chromosome) to sanity-check before running the whole genome.
- `--threads N` parallelizes across gene models; pick N = cores you have.
- The TSV is large for whole transcriptomes; `-o out.tsv.gz` compresses on the fly.
- `txtools makeDT --help` lists every option.
