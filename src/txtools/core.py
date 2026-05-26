"""Genomic -> transcriptomic read conversion (``tx_reads``).

This is the heart of txtools and one of the two performance bottlenecks. For
each gene model we:

1. fetch the alignments overlapping the gene's genomic span (BAM index),
2. pair mates (paired-end) and keep only reads/pairs fully contained in the
   gene's exons, on the matching strand,
3. validate that spliced reads' junctions coincide with annotated exon
   boundaries,
4. map each fragment's outer boundaries to transcriptomic coordinates, and
5. stitch the mate sequences into a single transcript-oriented sequence
   (``.`` fills the read1/read2 insert gap; overlaps are reconciled with an
   IUPAC consensus).

Boundary mapping uses the fragment's genomic extent together with the gene
strand rather than the read1/read2 mate roles. For ``strandMode = 1`` this is
numerically identical to the R package; for ``strandMode = 2`` it avoids the
inverted-coordinate behaviour of the original mate-role formula.
"""
from __future__ import annotations

import multiprocessing as mp
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np
import pysam

from .io import BamSource, GeneAnnotation, GeneModel
from .iupac import consensus_two

# CIGAR op codes (pysam): 0=M 1=I 2=D 3=N 4=S 5=H 6=P 7== 8=X
_M_OPS = (0, 7, 8)
_COMPLEMENT = str.maketrans("ACGTNacgtn-.", "TGCANtgcan-.")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _ref_len_no_N(cigartuples) -> int:
    """Aligned reference span with introns removed: sum of M/=/X and D lengths."""
    return sum(length for op, length in (cigartuples or ()) if op in (0, 2, 7, 8))


def _project_and_blocks(read: pysam.AlignedSegment):
    """Project a read's sequence into reference space and list its M-blocks.

    Returns ``(ref_seq, blocks)`` where ``ref_seq`` is the query laid onto the
    reference with introns (``N``) removed, deletions emitted as ``-`` and
    insertions/soft-clips dropped (equivalent to Biostrings
    ``sequenceLayer(to="reference-N-regions-removed")``); and ``blocks`` is a
    list of ``(start1, end1)`` 1-based inclusive reference ranges of each
    contiguous aligned (M/=/X) segment.
    """
    qseq = read.query_sequence or ""
    out = []
    blocks = []
    qpos = 0
    refpos = read.reference_start  # 0-based
    blk_start = None
    blk_end = None
    for op, length in read.cigartuples or ():
        if op in _M_OPS:
            out.append(qseq[qpos:qpos + length])
            if blk_start is None:
                blk_start = refpos + 1  # to 1-based
            blk_end = refpos + length    # 1-based inclusive
            qpos += length
            refpos += length
        elif op == 2:  # D: deletion -> gap, splits M-block
            out.append("-" * length)
            refpos += length
            # A deletion interrupts the M-run for junction-boundary purposes.
            if blk_start is not None:
                blocks.append((blk_start, blk_end))
                blk_start = None
        elif op == 3:  # N: intron, removed; splits block
            if blk_start is not None:
                blocks.append((blk_start, blk_end))
                blk_start = None
            refpos += length
        elif op == 1 or op == 4:  # I or S: consume query only
            qpos += length
        # H (5), P (6): consume nothing
    if blk_start is not None:
        blocks.append((blk_start, blk_end))
    return "".join(out), blocks


@dataclass
class GeneReads:
    """Transcriptomic reads assigned to one gene (analog of a GRanges entry)."""

    gene: GeneModel
    tx_start: np.ndarray          # 1-based tx coordinate of each fragment 5' end
    tx_end: np.ndarray            # 1-based tx coordinate of each fragment 3' end
    seqs: Optional[List[str]]     # merged transcript-oriented sequence, or None

    def __len__(self) -> int:
        return len(self.tx_start)


class TxReads:
    """Per-gene transcriptomic reads (analog of the R GRangesList)."""

    def __init__(self, by_gene: Dict[str, GeneReads], with_seq: bool):
        self.by_gene = by_gene
        self.with_seq = with_seq

    def __len__(self) -> int:
        return len(self.by_gene)

    def __iter__(self):
        return iter(self.by_gene.values())

    def __getitem__(self, name: str) -> GeneReads:
        return self.by_gene[name]

    @property
    def gene_names(self) -> List[str]:
        return list(self.by_gene.keys())

    @property
    def n_reads(self) -> int:
        return int(sum(len(g) for g in self.by_gene.values()))


def _passes_flags(read: pysam.AlignedSegment, load_secondary) -> bool:
    if read.is_unmapped or read.is_duplicate or read.is_qcfail or read.is_supplementary:
        return False
    if load_secondary is False and read.is_secondary:
        return False
    if load_secondary is True and not read.is_secondary:
        return False
    return True


def _junctions_ok(blocks, exon_starts: set, exon_ends: set) -> bool:
    """True if every M-block boundary coincides with an annotated exon edge.

    Single-block reads (no introns) always pass. Matches the intent of
    txtools' junction check, applied per-read (the original R code dropped all
    spliced reads in a gene if any one failed; here only the offending read is
    dropped).
    """
    if len(blocks) <= 1:
        return True
    for start1, end1 in blocks:
        if not (end1 in exon_ends or start1 in exon_starts):
            return False
    return True


def _merge_pair(seq5: str, seq3: Optional[str], width: int) -> str:
    """Stitch oriented 5' and 3' mate sequences into a length-``width`` string."""
    if seq5 and len(seq5) > width:
        seq5 = seq5[:width]
    if seq3 and len(seq3) > width:
        seq3 = seq3[-width:]
    merged = ["."] * width
    filled = bytearray(width)
    for i, ch in enumerate(seq5):
        merged[i] = ch
        filled[i] = 1
    if seq3:
        start3 = width - len(seq3)
        for j, ch in enumerate(seq3):
            pos = start3 + j
            if pos < 0:
                continue
            if filled[pos]:
                merged[pos] = consensus_two(merged[pos], ch)
            else:
                merged[pos] = ch
                filled[pos] = 1
    return "".join(merged)


def _process_gene_paired(gene: GeneModel, reads, with_seq: bool, strand_mode: int,
                         ignore_strand: bool):
    pos_to_tx = gene.pos_to_tx
    exon_starts = gene.exon_boundary_starts
    exon_ends = gene.exon_boundary_ends
    plus = gene.strand == "+"

    # Pair mates by query name within this region.
    mates: Dict[str, List[pysam.AlignedSegment]] = {}
    for r in reads:
        mates.setdefault(r.query_name, []).append(r)

    tx_starts, tx_ends, seqs = [], [], ([] if with_seq else None)
    for name, group in mates.items():
        r1 = r2 = None
        for r in group:
            if r.is_read1 and r1 is None:
                r1 = r
            elif r.is_read2 and r2 is None:
                r2 = r
        if r1 is None or r2 is None:
            continue

        # Pair strand (strandMode): mate whose mapped strand defines the pair.
        if strand_mode == 1:
            pair_minus = r1.is_reverse
        else:  # strand_mode == 2
            pair_minus = r2.is_reverse
        pair_strand = "-" if pair_minus else "+"
        if not ignore_strand and pair_strand != gene.strand:
            continue

        # 1-based outer boundaries of each mate.
        a_s, a_e = r1.reference_start + 1, r1.reference_end
        b_s, b_e = r2.reference_start + 1, r2.reference_end

        # All four boundaries must lie within exonic positions.
        if (a_s not in pos_to_tx or a_e not in pos_to_tx
                or b_s not in pos_to_tx or b_e not in pos_to_tx):
            continue

        proj1 = blocks1 = proj2 = blocks2 = None
        if r1.cigartuples and any(op == 3 for op, _ in r1.cigartuples):
            proj1, blocks1 = _project_and_blocks(r1)
            if not _junctions_ok(blocks1, exon_starts, exon_ends):
                continue
        if r2.cigartuples and any(op == 3 for op, _ in r2.cigartuples):
            proj2, blocks2 = _project_and_blocks(r2)
            if not _junctions_ok(blocks2, exon_starts, exon_ends):
                continue

        # Fragment genomic extent -> transcriptomic boundaries via strand.
        left = min(a_s, b_s)
        right = max(a_e, b_e)
        if plus:
            tx_s, tx_e = pos_to_tx[left], pos_to_tx[right]
        else:
            tx_s, tx_e = pos_to_tx[right], pos_to_tx[left]
        if tx_e < tx_s - 1:
            continue
        tx_starts.append(tx_s)
        tx_ends.append(tx_e)

        if with_seq:
            # Identify the 5' mate (defines transcript start) by position.
            if plus:
                five, three = (r1, r2) if a_s <= b_s else (r2, r1)
            else:
                five, three = (r1, r2) if a_e >= b_e else (r2, r1)
            if five is r1 and proj1 is not None:
                s5 = proj1
            elif five is r2 and proj2 is not None:
                s5 = proj2
            else:
                s5, _ = _project_and_blocks(five)
            if three is r1 and proj1 is not None:
                s3 = proj1
            elif three is r2 and proj2 is not None:
                s3 = proj2
            else:
                s3, _ = _project_and_blocks(three)
            if not plus:
                s5, s3 = _revcomp(s5), _revcomp(s3)
            width = tx_e - tx_s + 1
            seqs.append(_merge_pair(s5, s3, width))

    return _finalize(gene, tx_starts, tx_ends, seqs, with_seq)


def _process_gene_single(gene: GeneModel, reads, with_seq: bool, strand_mode: int,
                         ignore_strand: bool):
    pos_to_tx = gene.pos_to_tx
    exon_starts = gene.exon_boundary_starts
    exon_ends = gene.exon_boundary_ends
    plus = gene.strand == "+"

    tx_starts, tx_ends, seqs = [], [], ([] if with_seq else None)
    for r in reads:
        read_minus = r.is_reverse
        if strand_mode == 2:
            read_minus = not read_minus
        read_strand = "-" if read_minus else "+"
        if not ignore_strand and read_strand != gene.strand:
            continue

        s1, e1 = r.reference_start + 1, r.reference_end
        if s1 not in pos_to_tx or e1 not in pos_to_tx:
            continue

        proj = blocks = None
        if r.cigartuples and any(op == 3 for op, _ in r.cigartuples):
            proj, blocks = _project_and_blocks(r)
            if not _junctions_ok(blocks, exon_starts, exon_ends):
                continue

        if plus:
            tx_s, tx_e = pos_to_tx[s1], pos_to_tx[e1]
        else:
            tx_s, tx_e = pos_to_tx[e1], pos_to_tx[s1]
        if tx_e < tx_s - 1:
            continue
        # Drop reads whose transcriptomic footprint length differs from their
        # aligned reference length (introns removed): these skip an annotated
        # exon and don't fit this transcript model. Matches txtools' single-end
        # `width(tReads) == cigarWidth(N.regions.removed)` filter.
        if _ref_len_no_N(r.cigartuples) != (tx_e - tx_s + 1):
            continue
        tx_starts.append(tx_s)
        tx_ends.append(tx_e)

        if with_seq:
            if proj is None:
                proj, _ = _project_and_blocks(r)
            if not plus:
                proj = _revcomp(proj)
            width = tx_e - tx_s + 1
            # Drop reads whose projected length doesn't fit the tx footprint
            # (matches the single-end length check in txtools).
            seqs.append(_merge_pair(proj, None, width))

    return _finalize(gene, tx_starts, tx_ends, seqs, with_seq)


def _finalize(gene, tx_starts, tx_ends, seqs, with_seq) -> GeneReads:
    return GeneReads(
        gene=gene,
        tx_start=np.asarray(tx_starts, dtype=np.int64),
        tx_end=np.asarray(tx_ends, dtype=np.int64),
        seqs=seqs if with_seq else None,
    )


# --------------------------------------------------------------------------- #
# Worker plumbing for multiprocessing (each process opens its own BAM handle)
# --------------------------------------------------------------------------- #
_WORKER = {}


def _worker_init(bam_path: str):
    _WORKER["bam"] = pysam.AlignmentFile(bam_path, "rb")


def _worker_process(args):
    gene, paired, with_seq, strand_mode, ignore_strand, load_secondary = args
    bam = _WORKER["bam"]
    reads = [
        r for r in bam.fetch(gene.chrom, gene.start - 1, gene.end)
        if _passes_flags(r, load_secondary)
    ]
    if paired:
        gr = _process_gene_paired(gene, reads, with_seq, strand_mode, ignore_strand)
    else:
        gr = _process_gene_single(gene, reads, with_seq, strand_mode, ignore_strand)
    return gr


def tx_reads(
    reads: BamSource,
    gene_annot: GeneAnnotation,
    min_reads: int = 50,
    with_seq: bool = False,
    ignore_strand: bool = False,
    n_cores: int = 1,
    verbose: bool = True,
) -> TxReads:
    """Assign alignments to gene models and convert to transcriptomic space.

    Parameters mirror the R ``tx_reads``. ``reads`` is a :class:`BamSource`
    from :func:`txtools.io.tx_load_bam`. Genes retaining at least ``min_reads``
    assigned fragments are kept.
    """
    if with_seq and not reads.load_seq:
        raise ValueError(
            "with_seq=True requires the BAM to be loaded with load_seq=True "
            "(tx_load_bam(..., load_seq=True))."
        )

    # Restrict to genes on chromosomes present in the BAM.
    with reads.open() as bam:
        bam_chroms = set(bam.references)
    genes = [g for g in gene_annot if g.chrom in bam_chroms]
    if not genes:
        raise ValueError("No gene-model chromosomes are present in the BAM.")
    if verbose:
        print(f"Processing {len(genes)} gene models "
              f"({'paired' if reads.paired_end else 'single'}-end)...")

    tasks = [
        (g, reads.paired_end, with_seq, reads.strand_mode, ignore_strand,
         reads.load_secondary)
        for g in genes
    ]

    results: List[GeneReads] = []
    if n_cores and n_cores > 1:
        ctx = mp.get_context("fork") if hasattr(mp, "get_context") else mp
        with ctx.Pool(n_cores, initializer=_worker_init,
                      initargs=(reads.path,)) as pool:
            for gr in pool.imap_unordered(_worker_process, tasks, chunksize=16):
                results.append(gr)
    else:
        _worker_init(reads.path)
        for t in tasks:
            results.append(_worker_process(t))
        _WORKER.pop("bam", None)

    by_gene = {gr.gene.name: gr for gr in results if len(gr) >= min_reads}
    out = TxReads(by_gene, with_seq=with_seq)
    if verbose:
        print(f"Output: {out.n_reads} fragments in {len(out)} gene models.")
    return out
