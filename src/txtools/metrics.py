"""Derived per-position metrics (``tx_add_*``).

Each function takes a txDT (pandas DataFrame as produced by the ``tx_makeDT_*``
functions) and returns a new DataFrame with one extra column, mirroring the
txtools R API. Gene-wise operations (ratio shifts, rolling means, relative
position, motif search) respect gene boundaries; within each gene the rows are
contiguous and ordered by ``txcoor`` as produced by :mod:`txtools.tables`.
"""
from __future__ import annotations

import re
from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd

from .io import GeneAnnotation, Genome
from .iupac import IUPAC_CODE_2NUCS
from .tables import _ref_seq  # reuse strand-aware reference extraction

# Columns counted as "real" nucleotide reads (excludes N and the insert '.').
_TOTAL_NUCS = [n for n in IUPAC_CODE_2NUCS if n not in (".", "N")]
_MISINC_NUCS = ["A", "C", "G", "T", "-"]

_IUPAC_REGEX = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGTN]",
}


def _require(df: pd.DataFrame, cols: Iterable[str], who: str):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{who} requires column(s) {missing} in the txDT.")


def _drop_col(df: pd.DataFrame, col: str) -> pd.DataFrame:
    return df.drop(columns=[col]) if col in df.columns else df


def _gene_groups(df: pd.DataFrame):
    """Group by gene preserving table order (observed categories only)."""
    return df.groupby("gene", sort=False, observed=True)


def _insert_after(df: pd.DataFrame, col: str, after: str) -> pd.DataFrame:
    """Move ``col`` to immediately follow ``after`` (matches R's .after)."""
    cols = list(df.columns)
    cols.remove(col)
    idx = cols.index(after) + 1
    cols.insert(idx, col)
    return df[cols]


# --------------------------------------------------------------------------- #
# Read-start / read-end ratios
# --------------------------------------------------------------------------- #
def _ratio(df, num_col, min_cov):
    r = df[num_col] / df["cov"]
    r = r.where(df["cov"] >= min_cov, other=np.nan)
    return r.to_numpy(dtype=float)


def tx_add_startRatio(df: pd.DataFrame, min_cov: int = 50) -> pd.DataFrame:
    _require(df, ["start_5p", "cov"], "tx_add_startRatio")
    df = _drop_col(df, "startRatio").copy()
    df["startRatio"] = _ratio(df, "start_5p", min_cov)
    return df


def tx_add_endRatio(df: pd.DataFrame, min_cov: int = 50) -> pd.DataFrame:
    _require(df, ["end_3p", "cov"], "tx_add_endRatio")
    df = _drop_col(df, "endRatio").copy()
    df["endRatio"] = _ratio(df, "end_3p", min_cov)
    return df


def _shifted_ratio(df, num_col, min_cov, new_col, downstream):
    """Gene-wise ratio shifted 1 bp downstream (DS) or upstream (US).

    DS: each value takes the next position's value (last per gene -> NaN).
    US: each value takes the previous position's value (first per gene -> NaN).
    """
    df = _drop_col(df, new_col).copy()
    base = pd.Series(_ratio(df, num_col, min_cov), index=df.index)
    # shift(-1) brings the downstream value to the current row; shift(1) the upstream.
    shifted = _gene_groups(df.assign(_b=base))["_b"].shift(-1 if downstream else 1)
    df[new_col] = shifted.to_numpy(dtype=float)
    return df


def tx_add_startRatio1bpDS(df, min_cov=50):
    _require(df, ["start_5p", "cov"], "tx_add_startRatio1bpDS")
    return _shifted_ratio(df, "start_5p", min_cov, "startRatio1bpDS", downstream=True)


def tx_add_startRatio1bpUS(df, min_cov=50):
    _require(df, ["start_5p", "cov"], "tx_add_startRatio1bpUS")
    return _shifted_ratio(df, "start_5p", min_cov, "startRatio1bpUS", downstream=False)


def tx_add_endRatio1bpDS(df, min_cov=50):
    _require(df, ["end_3p", "cov"], "tx_add_endRatio1bpDS")
    return _shifted_ratio(df, "end_3p", min_cov, "endRatio1bpDS", downstream=True)


def tx_add_endRatio1bpUS(df, min_cov=50):
    _require(df, ["end_3p", "cov"], "tx_add_endRatio1bpUS")
    return _shifted_ratio(df, "end_3p", min_cov, "endRatio1bpUS", downstream=False)


# --------------------------------------------------------------------------- #
# Nucleotide totals / misincorporation
# --------------------------------------------------------------------------- #
def tx_add_nucTotal(df: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in _TOTAL_NUCS if c in df.columns]
    if not cols:
        raise ValueError("tx_add_nucTotal requires nucleotide-frequency columns.")
    df = _drop_col(df, "nucTotal").copy()
    df["nucTotal"] = df[cols].sum(axis=1).to_numpy()
    return df


def tx_add_misincCount(df: pd.DataFrame) -> pd.DataFrame:
    _require(df, ["refSeq"], "tx_add_misincCount")
    df = _drop_col(df, "misincCount").copy()
    present = [c for c in _MISINC_NUCS if c in df.columns]
    sub = df[present].to_numpy()
    ref = df["refSeq"].to_numpy().astype("U1")
    out = np.full(len(df), np.nan)
    for base in np.unique(ref):
        sel_cols = [i for i, c in enumerate(present) if c != base]
        mask = ref == base
        out[mask] = sub[np.ix_(mask, sel_cols)].sum(axis=1)
    df["misincCount"] = out
    return df


def tx_add_misincRate(df: pd.DataFrame, min_nuc_reads: int = 20,
                      add_counts: bool = False) -> pd.DataFrame:
    _require(df, ["refSeq"], "tx_add_misincRate")
    df = _drop_col(df, "misincRate").copy()
    tmp = tx_add_nucTotal(tx_add_misincCount(df))
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.round(tmp["misincCount"].to_numpy() / tmp["nucTotal"].to_numpy(), 6)
    rate[tmp["nucTotal"].to_numpy() < min_nuc_reads] = np.nan
    if add_counts:
        df = tmp
    df["misincRate"] = rate
    return df


def tx_add_misincRateNucSpec(df: pd.DataFrame, ref_nuc: str, mis_nuc: str,
                             min_nuc_reads: int = 20) -> pd.DataFrame:
    _require(df, ["refSeq", ref_nuc, mis_nuc], "tx_add_misincRateNucSpec")
    col = f"MR_{ref_nuc}to{mis_nuc}"
    df = _drop_col(df, col).copy()
    mis = df[mis_nuc].to_numpy(dtype=float)
    denom = mis + df[ref_nuc].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.round(mis / denom, 6)
    rate[denom < min_nuc_reads] = np.nan
    rate[df["refSeq"].to_numpy().astype("U1") != ref_nuc] = np.nan
    df[col] = rate
    return df


# --------------------------------------------------------------------------- #
# Annotation columns
# --------------------------------------------------------------------------- #
def tx_add_geneRegion(df: pd.DataFrame, gene_annot: GeneAnnotation) -> pd.DataFrame:
    _require(df, ["gencoor", "strand", "gene"], "tx_add_geneRegion")
    df = _drop_col(df, "geneRegion").copy()
    region = np.empty(len(df), dtype=object)
    gencoor = df["gencoor"].to_numpy()
    genes = df["gene"].to_numpy().astype(str)
    for name in pd.unique(genes):
        g = gene_annot[name]
        mask = genes == name
        if g.thick_end < g.thick_start:  # non-coding (zero-width CDS)
            region[mask] = "non-coding"
            continue
        gc = gencoor[mask]
        lab = np.empty(gc.shape, dtype=object)
        in_cds = (gc >= g.thick_start) & (gc <= g.thick_end)
        if g.strand == "+":
            lab[gc < g.thick_start] = "5'UTR"
            lab[gc > g.thick_end] = "3'UTR"
        else:
            lab[gc > g.thick_end] = "5'UTR"
            lab[gc < g.thick_start] = "3'UTR"
        lab[in_cds] = "CDS"
        region[mask] = lab
    df["geneRegion"] = pd.Categorical(
        region, categories=["5'UTR", "CDS", "3'UTR", "non-coding"])
    return df


def _motif_regex(motif: str) -> re.Pattern:
    body = "".join(_IUPAC_REGEX[ch.upper()] for ch in motif)
    # Lookahead enables overlapping matches, matching Biostrings matchPattern.
    return re.compile(f"(?=({body}))")


def tx_add_motifPresence(df: pd.DataFrame, motif: str, nuc_positions="all",
                         motif_col_name: str = "auto",
                         mask_N: bool = True) -> pd.DataFrame:
    _require(df, ["refSeq", "gene"], "tx_add_motifPresence")
    if motif_col_name == "auto":
        motif_col_name = f"{motif}_motif_{nuc_positions}"
    df = _drop_col(df, motif_col_name).copy()

    mid = None
    if nuc_positions == "center":
        if len(motif) % 2 == 0:
            raise ValueError("Even-length motif has no single center position.")
        mid = (len(motif) + 1) // 2
    elif isinstance(nuc_positions, int):
        if nuc_positions > len(motif):
            raise ValueError("nuc_positions exceeds motif length.")
    elif nuc_positions != "all":
        raise ValueError("nuc_positions must be 'all', 'center', or an int.")

    rx = _motif_regex(motif)
    flags = np.zeros(len(df), dtype=bool)
    ref = df["refSeq"].to_numpy().astype("U1")
    genes = df["gene"].to_numpy().astype(str)

    start = 0
    n = len(df)
    while start < n:
        gname = genes[start]
        end = start
        while end < n and genes[end] == gname:
            end += 1
        seq = "".join(ref[start:end])
        if mask_N:
            seq = seq.replace("N", ".")
        for m in rx.finditer(seq):
            s = m.start()
            mlen = len(m.group(1))
            if nuc_positions == "all":
                flags[start + s: start + s + mlen] = True
            elif nuc_positions == "center":
                flags[start + s + mid - 1] = True
            else:  # int
                flags[start + s + nuc_positions - 1] = True
        start = end

    df[motif_col_name] = flags
    return df


def tx_add_rollingMean(df: pd.DataFrame, col_name: str, win_size: int,
                       new_col_name: Optional[str] = None, align: str = "center",
                       min_cov: int = 20) -> pd.DataFrame:
    _require(df, [col_name], "tx_add_rollingMean")
    if new_col_name is None:
        new_col_name = f"rollMean_{col_name}_{win_size}"
    df = _drop_col(df, new_col_name).copy()

    def roll(s: pd.Series) -> pd.Series:
        if align == "center":
            return s.rolling(win_size, center=True, min_periods=win_size).mean()
        if align == "right":
            return s.rolling(win_size, min_periods=win_size).mean()
        if align == "left":
            return s.rolling(win_size, min_periods=win_size).mean().shift(-(win_size - 1))
        raise ValueError("align must be 'left', 'center' or 'right'.")

    rolled = _gene_groups(df)[col_name].transform(roll)
    out = rolled.to_numpy(dtype=float)
    if "cov" in df.columns:
        out[df["cov"].to_numpy() <= min_cov] = np.nan
    df[new_col_name] = out
    return df


def tx_add_relTxPos(df: pd.DataFrame, round_dig: int = 3) -> pd.DataFrame:
    _require(df, ["txcoor"], "tx_add_relTxPos")
    df = _drop_col(df, "relTxPos").copy()
    maxc = _gene_groups(df)["txcoor"].transform("max")
    df["relTxPos"] = np.round(df["txcoor"].to_numpy() / maxc.to_numpy(), round_dig)
    return df


def tx_add_pos(df: pd.DataFrame, sep: str = ":", check_uniq: bool = True) -> pd.DataFrame:
    _require(df, ["gene", "txcoor"], "tx_add_pos")
    df = _drop_col(df, "pos").copy()
    pos = df["gene"].astype(str) + sep + df["txcoor"].astype(str)
    if check_uniq and pos.duplicated().any():
        raise ValueError("gene + txcoor combinations are not unique.")
    df["pos"] = pos.to_numpy()
    return _insert_after(df, "pos", "txcoor")


def tx_add_refSeqDT(df: pd.DataFrame, genome: Genome,
                    gene_annot: GeneAnnotation) -> pd.DataFrame:
    """Add (or replace) the strand-aware reference-sequence column."""
    _require(df, ["gene", "txcoor"], "tx_add_refSeqDT")
    df = _drop_col(df, "refSeq").copy()
    genes = df["gene"].to_numpy().astype(str)
    ref = np.empty(len(df), dtype="U1")
    start, n = 0, len(df)
    while start < n:
        gname = genes[start]
        end = start
        while end < n and genes[end] == gname:
            end += 1
        ref[start:end] = _ref_seq(gene_annot[gname], genome)
        start = end
    df["refSeq"] = ref
    return _insert_after(df, "refSeq", "txcoor")


def tx_add_siteAnnotation(df: pd.DataFrame, sites: Sequence, col_name: str) -> pd.DataFrame:
    """Mark 1-nt sites given as ``(chrom, gencoor, strand)`` tuples.

    Replacement for the GRanges-based ``tx_add_siteAnnotation``; accepts an
    iterable of ``(chrom, gencoor, strand)`` (1-based genomic coordinate).
    """
    _require(df, ["chr", "gencoor", "strand"], "tx_add_siteAnnotation")
    df = _drop_col(df, col_name).copy()
    site_set = {(str(c), int(p), str(s)) for c, p, s in sites}
    key = list(zip(df["chr"].astype(str), df["gencoor"].astype(int),
                   df["strand"].astype(str)))
    df[col_name] = np.fromiter((k in site_set for k in key), dtype=bool, count=len(df))
    return df
