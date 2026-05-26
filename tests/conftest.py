"""Shared test fixtures: a self-contained toy genome / annotation / BAM.

Builds a small dataset on disk exercising:
* a ``+`` strand two-exon gene with a spliced read and a known mismatch;
* a ``-`` strand single-exon gene with a perfect-match pair.
The transcript sequences are recorded so tests can assert exact values.
"""
import os

import numpy as np
import pysam
import pytest


def _revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _make_read(name, is_read1, reverse, ref_start0, qseq, cigar, mate_start0):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = qseq
    a.is_paired = True
    a.is_proper_pair = True
    a.is_read1 = is_read1
    a.is_read2 = not is_read1
    a.is_reverse = reverse
    a.mate_is_reverse = not reverse
    a.reference_id = 0
    a.reference_start = ref_start0
    a.next_reference_id = 0
    a.next_reference_start = mate_start0
    a.mapping_quality = 60
    a.cigartuples = cigar
    return a


@pytest.fixture(scope="session")
def toy(tmp_path_factory):
    d = tmp_path_factory.mktemp("txtools_toy")
    chrseq = "".join(np.random.RandomState(7).choice(list("ACGT"), 300))

    fa = os.path.join(d, "toy.fa")
    with open(fa, "w") as fh:
        fh.write(">chr1 test genome\n")
        for i in range(0, len(chrseq), 60):
            fh.write(chrseq[i:i + 60] + "\n")

    bed = os.path.join(d, "toy.bed")
    with open(bed, "w") as fh:
        # gene1 '+': two exons, genomic 1-based 11..40 and 61..90  (tx len 60)
        fh.write("chr1\t10\t90\tgene1\t0\t+\t10\t90\t0\t2\t30,30\t0,50\n")
        # gene2 '-': single exon, genomic 1-based 101..150          (tx len 50)
        fh.write("chr1\t100\t150\tgene2\t0\t-\t100\t150\t0\t1\t50\t0\n")

    tx_gene1 = chrseq[10:40] + chrseq[60:90]
    tx_gene2 = _revcomp(chrseq[100:150])

    reads = []
    # gene1 pair A: read1 fwd 11..30 (tx1..20), read2 rev 71..90 (tx41..60)
    reads.append(_make_read("A", True, False, 10, chrseq[10:30], [(0, 20)], 70))
    reads.append(_make_read("A", False, True, 70, chrseq[70:90], [(0, 20)], 10))
    # gene1 pair B: spliced read1 31..40 |N| 61..70, read2 rev 76..90
    seqB1 = chrseq[30:40] + chrseq[60:70]
    reads.append(_make_read("B", True, False, 30, seqB1, [(0, 10), (3, 20), (0, 10)], 75))
    reads.append(_make_read("B", False, True, 75, chrseq[75:90], [(0, 15)], 30))
    # gene1 pair C: like A but mismatch at tx5 (genomic 15) C/A -> introduce change
    seqC1 = list(chrseq[10:30])
    seqC1[4] = "A" if seqC1[4] != "A" else "C"
    reads.append(_make_read("C", True, False, 10, "".join(seqC1), [(0, 20)], 70))
    reads.append(_make_read("C", False, True, 70, chrseq[70:90], [(0, 20)], 10))

    # gene2 '-' pair D: read1 rev 131..150 (tx1..20), read2 fwd 101..120 (tx31..50)
    reads.append(_make_read("D", True, True, 130, chrseq[130:150], [(0, 20)], 100))
    reads.append(_make_read("D", False, False, 100, chrseq[100:120], [(0, 20)], 130))
    # gene2 pair E: same geometry, perfect match (gives min_reads=2)
    reads.append(_make_read("E", True, True, 130, chrseq[130:150], [(0, 20)], 100))
    reads.append(_make_read("E", False, False, 100, chrseq[100:120], [(0, 20)], 130))

    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": len(chrseq)}]}
    bam = os.path.join(d, "toy.bam")
    reads.sort(key=lambda r: r.reference_start)
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.index(bam)

    return {
        "dir": str(d), "fasta": fa, "bed": bed, "bam": bam,
        "chrseq": chrseq, "tx_gene1": tx_gene1, "tx_gene2": tx_gene2,
    }
