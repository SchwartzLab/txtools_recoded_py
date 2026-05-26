"""Summarized nucleotide-resolution tables (``tx_makeDT_*``).

Builds the txDT: one row per transcriptomic position per gene, with the core
coordinate columns (``chr, gencoor, strand, gene, txcoor``) plus, depending on
the function, coverage metrics and/or a nucleotide-frequency pileup, and an
optional reference-sequence column.

The pileup is the second performance bottleneck, so coverage uses a
difference-array (O(reads + L)) and the nucleotide tally maps sequence
characters to matrix columns through a 256-entry byte table, keeping the inner
loops in NumPy/C rather than Python.
"""
from __future__ import annotations

from typing import List, Optional

import numpy as np
import pandas as pd

from .core import GeneReads, TxReads
from .io import GeneAnnotation, Genome, GeneModel
from .iupac import (
    IUPAC_CODE_2NUCS,
    SYMBOL_INDEX,
    matrix_columns,
    simplify_matrix,
)

CORE_COLS = ["chr", "gencoor", "strand", "gene", "txcoor"]

# 256-entry translation table: byte value of each character -> matrix column.
_N_COL = SYMBOL_INDEX["N"]
_COL_TABLE = bytes(SYMBOL_INDEX.get(chr(i), _N_COL) for i in range(256))
_N_SYMBOLS = len(IUPAC_CODE_2NUCS)


# --------------------------------------------------------------------------- #
# Per-gene metric computation
# --------------------------------------------------------------------------- #
def _coverage_arrays(gr: GeneReads, length: int):
    """Return (cov, start_5p, end_3p) integer arrays of length ``length``."""
    txs = gr.tx_start
    txe = gr.tx_end
    if len(txs) == 0:
        z = np.zeros(length, dtype=np.int64)
        return z, z.copy(), z.copy()
    diff = np.zeros(length + 1, dtype=np.int64)
    np.add.at(diff, txs - 1, 1)
    np.add.at(diff, txe, -1)
    cov = np.cumsum(diff[:length])
    start_5p = np.bincount(txs - 1, minlength=length)[:length].astype(np.int64)
    end_3p = np.bincount(txe - 1, minlength=length)[:length].astype(np.int64)
    return cov, start_5p, end_3p


def _nucfreq_matrix(gr: GeneReads, length: int) -> np.ndarray:
    """Return an ``L x len(IUPAC_CODE_2NUCS)`` integer pileup matrix."""
    mat = np.zeros((length, _N_SYMBOLS), dtype=np.int64)
    if gr.seqs is None:
        raise ValueError(
            "Nucleotide frequency requires sequences; run tx_reads(with_seq=True)."
        )
    txs = gr.tx_start
    for i, seq in enumerate(gr.seqs):
        s0 = int(txs[i]) - 1
        cols = np.frombuffer(seq.encode("latin1").translate(_COL_TABLE),
                             dtype=np.uint8).astype(np.intp)
        if cols.size == 0:
            continue
        positions = np.arange(s0, s0 + cols.size)
        np.add.at(mat, (positions, cols), 1)
    return mat


def _ref_seq(gene: GeneModel, genome: Genome) -> np.ndarray:
    """Reference base per transcriptomic position (strand-aware)."""
    parts = [genome.fetch(gene.chrom, s, e)
             for s, e in zip(gene.exon_starts, gene.exon_ends)]
    seq = "".join(parts)
    if gene.strand == "-":
        from .core import _revcomp
        seq = _revcomp(seq)
    return np.frombuffer(seq.encode("latin1"), dtype="S1").astype("U1")


def _coord_columns(gene: GeneModel):
    length = gene.length
    gencoor = gene.exon_positions.astype(np.int64)
    txcoor = np.arange(1, length + 1, dtype=np.int64)
    return length, gencoor, txcoor


# --------------------------------------------------------------------------- #
# Table assembly
# --------------------------------------------------------------------------- #
def _genes_in_order(x: TxReads, gene_annot: GeneAnnotation, full: bool) -> List[GeneModel]:
    if full:
        return list(gene_annot)
    present = set(x.gene_names)
    return [g for g in gene_annot if g.name in present]


def _empty_reads(gene: GeneModel, with_seq: bool) -> GeneReads:
    return GeneReads(gene=gene,
                     tx_start=np.empty(0, dtype=np.int64),
                     tx_end=np.empty(0, dtype=np.int64),
                     seqs=([] if with_seq else None))


def _assemble(x: TxReads, gene_annot: GeneAnnotation, genome: Optional[Genome],
              full: bool, want_cov: bool, want_nuc: bool, simplify: str) -> pd.DataFrame:
    genes = _genes_in_order(x, gene_annot, full)
    if not genes:
        raise ValueError("No genes to tabulate.")

    chr_parts, gencoor_parts, strand_parts, gene_parts, txcoor_parts = [], [], [], [], []
    ref_parts = [] if genome is not None else None
    cov_parts, sta_parts, end_parts = [], [], []
    nuc_cols = matrix_columns(simplify) if want_nuc else []
    nuc_parts = {c: [] for c in nuc_cols}

    for gene in genes:
        length, gencoor, txcoor = _coord_columns(gene)
        gr = x.by_gene.get(gene.name) or _empty_reads(gene, x.with_seq)

        chr_parts.append(np.full(length, gene.chrom, dtype=object))
        gencoor_parts.append(gencoor)
        strand_parts.append(np.full(length, gene.strand, dtype=object))
        gene_parts.append(np.full(length, gene.name, dtype=object))
        txcoor_parts.append(txcoor)
        if ref_parts is not None:
            ref_parts.append(_ref_seq(gene, genome))

        if want_cov:
            cov, sta, end = _coverage_arrays(gr, length)
            cov_parts.append(cov)
            sta_parts.append(sta)
            end_parts.append(end)

        if want_nuc:
            mat = simplify_matrix(_nucfreq_matrix(gr, length), simplify)
            for j, c in enumerate(nuc_cols):
                nuc_parts[c].append(mat[:, j])

    data = {
        "chr": np.concatenate(chr_parts),
        "gencoor": np.concatenate(gencoor_parts),
        "strand": np.concatenate(strand_parts),
        "gene": np.concatenate(gene_parts),
        "txcoor": np.concatenate(txcoor_parts),
    }
    if ref_parts is not None:
        data["refSeq"] = np.concatenate(ref_parts)
    if want_cov:
        data["cov"] = np.concatenate(cov_parts)
        data["start_5p"] = np.concatenate(sta_parts)
        data["end_3p"] = np.concatenate(end_parts)
    if want_nuc:
        for c in nuc_cols:
            data[c] = np.concatenate(nuc_parts[c])

    df = pd.DataFrame(data)
    df["chr"] = df["chr"].astype("category")
    df["strand"] = df["strand"].astype("category")
    df["gene"] = df["gene"].astype("category")
    return df


# --------------------------------------------------------------------------- #
# Public makeDT functions (names kept consistent with the R package)
# --------------------------------------------------------------------------- #
def tx_makeDT_coverage(x: TxReads, gene_annot: GeneAnnotation,
                       genome: Optional[Genome] = None,
                       full_dt: bool = False) -> pd.DataFrame:
    """Coverage table: cov (insert coverage), start_5p, end_3p per position."""
    return _assemble(x, gene_annot, genome, full_dt,
                     want_cov=True, want_nuc=False, simplify="splitForceInt")


def tx_makeDT_nucFreq(x: TxReads, gene_annot: GeneAnnotation,
                      genome: Optional[Genome] = None,
                      simplify_IUPAC: str = "splitForceInt",
                      full_dt: bool = False) -> pd.DataFrame:
    """Nucleotide-frequency table: per-base counts incl. deletions and inserts."""
    return _assemble(x, gene_annot, genome, full_dt,
                     want_cov=False, want_nuc=True, simplify=simplify_IUPAC)


def tx_makeDT_covNucFreq(x: TxReads, gene_annot: GeneAnnotation,
                         genome: Optional[Genome] = None,
                         simplify_IUPAC: str = "splitForceInt",
                         full_dt: bool = False) -> pd.DataFrame:
    """Combined coverage + nucleotide-frequency table."""
    return _assemble(x, gene_annot, genome, full_dt,
                     want_cov=True, want_nuc=True, simplify=simplify_IUPAC)
