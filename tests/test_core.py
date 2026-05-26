import numpy as np

import txtools


def _load(toy, **kw):
    ga = txtools.tx_load_bed(toy["bed"])
    genome = txtools.tx_load_genome(toy["fasta"])
    bam = txtools.tx_load_bam(toy["bam"], paired_end=True, load_seq=True, verbose=False)
    reads = txtools.tx_reads(bam, ga, min_reads=1, with_seq=True, verbose=False, **kw)
    return ga, genome, reads


def test_bed_exon_model(toy):
    ga = txtools.tx_load_bed(toy["bed"])
    g1 = ga["gene1"]
    assert g1.strand == "+"
    assert g1.length == 60
    # exon positions are 1-based genomic, 5'->3'
    assert g1.exon_positions[0] == 11
    assert g1.exon_positions[29] == 40
    assert g1.exon_positions[30] == 61
    assert g1.exon_positions[-1] == 90
    g2 = ga["gene2"]
    assert g2.strand == "-"
    assert g2.length == 50
    assert g2.exon_positions[0] == 150  # 5' end of a minus-strand gene
    assert g2.exon_positions[-1] == 101


def test_tx_reads_plus_strand_coords(toy):
    _, _, reads = _load(toy)
    gr = reads["gene1"]
    assert len(gr) == 3
    # pairs A and C span the whole transcript; B starts at the spliced junction tail
    assert sorted(gr.tx_start.tolist()) == [1, 1, 21]
    assert gr.tx_end.tolist() == [60, 60, 60]


def test_spliced_read_sequence(toy):
    _, _, reads = _load(toy)
    gr = reads["gene1"]
    # find the spliced fragment (tx_start == 21)
    i = int(np.where(gr.tx_start == 21)[0][0])
    seq = gr.seqs[i]
    assert len(seq) == 40                      # tx 21..60
    assert seq[20:25] == "....."               # read1/read2 insert gap (tx41..45)
    assert "." not in seq[:20]                  # exon-spanning portion is contiguous


def test_minus_strand_coords_and_refseq(toy):
    ga, genome, reads = _load(toy)
    gr = reads["gene2"]
    assert len(gr) == 2
    assert gr.tx_start.tolist() == [1, 1]
    assert gr.tx_end.tolist() == [50, 50]
    dt = txtools.tx_makeDT_covNucFreq(reads, ga, genome=genome)
    g2 = dt[dt.gene == "gene2"]
    assert "".join(g2.refSeq.tolist()) == toy["tx_gene2"]


def test_coverage_and_insert(toy):
    ga, genome, reads = _load(toy)
    dt = txtools.tx_makeDT_covNucFreq(reads, ga, genome=genome)
    g1 = dt[dt.gene == "gene1"].reset_index(drop=True)
    # tx1..20 covered by A and C only; tx21+ also by B
    assert g1.loc[0, "cov"] == 2
    assert g1.loc[20, "cov"] == 3
    # read-starts: A,C start at tx1; B starts at tx21
    assert g1.loc[0, "start_5p"] == 2
    assert g1.loc[20, "start_5p"] == 1
    # A and C have their read1/read2 insert gap at tx21..40 (read1 ends tx20,
    # read2 starts tx41); B carries real bases there. So '.'==2 in the interior.
    gap_AC = g1[g1.txcoor.between(22, 39)]
    assert (gap_AC["."] == 2).all()
    # B's own insert gap is tx41..45, where A and C carry real bases -> '.'==1.
    gap_B = g1[g1.txcoor.between(41, 45)]
    assert (gap_B["."] == 1).all()


def test_mismatch_pileup(toy):
    ga, genome, reads = _load(toy)
    dt = txtools.tx_makeDT_covNucFreq(reads, ga, genome=genome)
    row = dt[(dt.gene == "gene1") & (dt.txcoor == 5)].iloc[0]
    ref = row["refSeq"]
    # one read keeps the reference base, the other was mutated
    assert row[ref] == 1
    assert row.drop(labels=["chr", "gencoor", "strand", "gene", "txcoor", "refSeq"])[
        ["A", "C", "G", "T"]].sum() == 2


def test_ref_len_no_N():
    from txtools.core import _ref_len_no_N
    # ops: 0=M 1=I 2=D 3=N 4=S 7== 8=X ; reference span excludes I/S/H and N
    assert _ref_len_no_N([(0, 50)]) == 50                       # 50M
    assert _ref_len_no_N([(0, 50), (2, 2), (0, 49)]) == 101      # 50M2D49M (D counts)
    assert _ref_len_no_N([(0, 50), (3, 1000), (0, 51)]) == 101   # 50M1000N51M (N excluded)
    assert _ref_len_no_N([(4, 5), (0, 96)]) == 96                # 5S96M (soft-clip excluded)
    assert _ref_len_no_N([(0, 50), (1, 3), (0, 48)]) == 98       # 50M3I48M (insertion excluded)


def test_refseq_matches_transcript(toy):
    ga, genome, reads = _load(toy)
    dt = txtools.tx_makeDT_covNucFreq(reads, ga, genome=genome)
    g1 = dt[dt.gene == "gene1"]
    assert "".join(g1.refSeq.tolist()) == toy["tx_gene1"]
