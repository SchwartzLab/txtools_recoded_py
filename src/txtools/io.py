"""Loading input data: BAM alignments, BED gene models, FASTA genome.

These correspond to the ``tx_load_*`` family of the txtools R package. The
public entry points keep their original names (:func:`tx_load_bam`,
:func:`tx_load_bed`, :func:`tx_load_genome`) but return lightweight Python
objects instead of Bioconductor S4 objects.

Coordinate convention
----------------------
To stay numerically faithful to the R package (which is built on 1-based,
inclusive GRanges), all *genomic* coordinates exposed by these objects are
**1-based inclusive**. pysam's 0-based positions are converted at the
boundary. Transcriptomic coordinates (``txcoor``) are 1-based, 5'->3'.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
import pysam


# --------------------------------------------------------------------------- #
# Genome
# --------------------------------------------------------------------------- #
class Genome:
    """Random-access reference genome backed by a FASTA (``.fai`` indexed).

    Chromosome names are truncated to the first whitespace-delimited token,
    matching ``tx_load_genome``.
    """

    def __init__(self, fasta_path: str):
        self.path = fasta_path
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(fasta_path)
        fai = fasta_path + ".fai"
        if not os.path.exists(fai):
            pysam.faidx(fasta_path)  # creates the .fai next to the FASTA
        self._fa = pysam.FastaFile(fasta_path)
        # Map possibly-long FASTA names to their first token, and back.
        self._full_names = list(self._fa.references)
        self._short_to_full: Dict[str, str] = {}
        for full in self._full_names:
            short = full.split()[0]
            self._short_to_full.setdefault(short, full)

    @property
    def chromosomes(self) -> List[str]:
        return list(self._short_to_full.keys())

    def fetch(self, chrom: str, start1: int, end1: int) -> str:
        """Return uppercase sequence for 1-based inclusive ``[start1, end1]``."""
        full = self._short_to_full.get(chrom, chrom)
        return self._fa.fetch(full, start1 - 1, end1).upper()

    def close(self):
        self._fa.close()


def tx_load_genome(fasta_file: str) -> Genome:
    """Load a genome from a FASTA file (creates a ``.fai`` index if needed)."""
    return Genome(fasta_file)


# --------------------------------------------------------------------------- #
# Gene models (BED12 / BED6)
# --------------------------------------------------------------------------- #
@dataclass
class GeneModel:
    """A single transcript / gene model.

    All genomic coordinates are 1-based inclusive.
    """

    name: str
    chrom: str
    strand: str
    start: int            # 1-based genomic start of the gene span
    end: int              # 1-based genomic end of the gene span
    exon_starts: List[int]  # 1-based genomic starts of each exon (genome order)
    exon_ends: List[int]    # 1-based genomic ends of each exon (genome order)
    thick_start: int      # CDS start, 1-based (== start if non-coding)
    thick_end: int        # CDS end, 1-based (== start-1 / 0-width if non-coding)

    # Lazily computed transcriptomic mapping.
    _exon_pos: Optional[np.ndarray] = field(default=None, repr=False)
    _pos_to_tx: Optional[Dict[int, int]] = field(default=None, repr=False)

    @property
    def length(self) -> int:
        return int(sum(e - s + 1 for s, e in zip(self.exon_starts, self.exon_ends)))

    @property
    def exon_positions(self) -> np.ndarray:
        """1-based genomic positions of every exonic base, in 5'->3' order.

        This is the txtools ``iExon`` vector: index ``i`` (0-based) corresponds
        to transcriptomic coordinate ``i + 1``.
        """
        if self._exon_pos is None:
            parts = [
                np.arange(s, e + 1, dtype=np.int64)
                for s, e in zip(self.exon_starts, self.exon_ends)
            ]
            pos = np.concatenate(parts) if parts else np.empty(0, dtype=np.int64)
            if self.strand == "-":
                pos = pos[::-1].copy()
            self._exon_pos = pos
        return self._exon_pos

    @property
    def pos_to_tx(self) -> Dict[int, int]:
        """Map 1-based genomic position -> 1-based transcriptomic coordinate."""
        if self._pos_to_tx is None:
            self._pos_to_tx = {
                int(g): i + 1 for i, g in enumerate(self.exon_positions)
            }
        return self._pos_to_tx

    @property
    def exon_boundary_starts(self) -> set:
        return set(self.exon_starts)

    @property
    def exon_boundary_ends(self) -> set:
        return set(self.exon_ends)


class GeneAnnotation:
    """An ordered collection of :class:`GeneModel` keyed by gene name."""

    def __init__(self, genes: List[GeneModel]):
        self.genes = genes
        self._by_name = {g.name: g for g in genes}
        if len(self._by_name) != len(genes):
            dups = _duplicates([g.name for g in genes])
            raise ValueError(
                "Duplicated genes found in gene annotation: " + " ".join(dups)
            )

    def __len__(self) -> int:
        return len(self.genes)

    def __iter__(self):
        return iter(self.genes)

    def __getitem__(self, name: str) -> GeneModel:
        return self._by_name[name]

    @property
    def names(self) -> List[str]:
        return [g.name for g in self.genes]

    @property
    def chromosomes(self) -> List[str]:
        return sorted({g.chrom for g in self.genes})

    def subset(self, names) -> "GeneAnnotation":
        wanted = set(names)
        return GeneAnnotation([g for g in self.genes if g.name in wanted])


def _duplicates(items):
    seen, dups = set(), []
    for x in items:
        if x in seen and x not in dups:
            dups.append(x)
        seen.add(x)
    return dups


def tx_load_bed(bedfile: str) -> GeneAnnotation:
    """Load gene models from a BED-12 or BED-6 file.

    BED-12 preserves the exon (block) structure; BED-6 (or narrower) treats
    each record as a single-exon gene. Matches ``tx_load_bed`` including its
    rejection of duplicated gene names.
    """
    genes: List[GeneModel] = []
    with open(bedfile) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if len(f) < 3:
                raise ValueError(f"{bedfile}:{lineno}: fewer than 3 BED columns")
            chrom = f[0]
            chrom_start = int(f[1])           # 0-based
            chrom_end = int(f[2])             # exclusive
            name = f[3] if len(f) > 3 and f[3] != "" else f"{chrom}:{chrom_start}-{chrom_end}"
            strand = f[5] if len(f) > 5 and f[5] in "+-" else "+"
            start1 = chrom_start + 1          # 1-based inclusive
            end1 = chrom_end

            if len(f) >= 8:
                thick_start = int(f[6]) + 1
                thick_end = int(f[7])
            else:
                thick_start, thick_end = start1, start1 - 1  # no CDS

            if len(f) >= 12 and int(f[9]) > 0:
                block_count = int(f[9])
                sizes = [int(x) for x in f[10].rstrip(",").split(",")][:block_count]
                offsets = [int(x) for x in f[11].rstrip(",").split(",")][:block_count]
                exon_starts, exon_ends = [], []
                for off, sz in zip(offsets, sizes):
                    es = chrom_start + off + 1          # 1-based
                    exon_starts.append(es)
                    exon_ends.append(es + sz - 1)
            else:
                exon_starts = [start1]
                exon_ends = [end1]

            genes.append(
                GeneModel(
                    name=name, chrom=chrom, strand=strand,
                    start=start1, end=end1,
                    exon_starts=exon_starts, exon_ends=exon_ends,
                    thick_start=thick_start, thick_end=thick_end,
                )
            )
    return GeneAnnotation(genes)


# --------------------------------------------------------------------------- #
# BAM alignments
# --------------------------------------------------------------------------- #
@dataclass
class BamSource:
    """Handle to an indexed BAM, plus how its reads should be interpreted.

    Unlike ``tx_load_bam`` in R (which reads every alignment into memory), this
    keeps a reference to the file and its parameters; :func:`tx_reads` then
    fetches reads per gene region using the BAM index. This is far more memory
    efficient for whole-transcriptome runs while exposing the same call.
    """

    path: str
    paired_end: bool
    load_seq: bool = False
    strand_mode: int = 1
    load_secondary: Optional[bool] = False  # False=primary only, True=secondary only, None=both

    def open(self) -> pysam.AlignmentFile:
        return pysam.AlignmentFile(self.path, "rb")


def tx_load_bam(
    file: str,
    paired_end: bool,
    load_seq: bool = False,
    strand_mode: int = 1,
    load_secondary: Optional[bool] = False,
    verbose: bool = True,
) -> BamSource:
    """Prepare a BAM file for transcriptomic processing.

    Parameters mirror the R ``tx_load_bam``. The BAM must be coordinate-sorted
    and indexed; a ``.bai`` is created automatically if missing.
    """
    if not os.path.exists(file):
        raise FileNotFoundError(file)
    if strand_mode not in (1, 2):
        raise ValueError("strand_mode must be 1 or 2 (mode 0/'*' is unsupported)")

    with pysam.AlignmentFile(file, "rb") as bam:
        so = bam.header.to_dict().get("HD", {}).get("SO")
        if so != "coordinate":
            raise ValueError(
                f"BAM {file} is not coordinate-sorted (SO={so!r}). "
                "Sort it first: `samtools sort`."
            )
    if not (os.path.exists(file + ".bai") or os.path.exists(file[:-1] + "i")
            or os.path.exists(os.path.splitext(file)[0] + ".bai")):
        if verbose:
            print(f"Indexing BAM (no .bai found): {file}")
        pysam.index(file)

    return BamSource(
        path=file,
        paired_end=paired_end,
        load_seq=load_seq,
        strand_mode=strand_mode,
        load_secondary=load_secondary,
    )
