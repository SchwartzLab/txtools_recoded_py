import numpy as np

import txtools
from txtools.iupac import simplify_force_int, IUPAC_CODE_2NUCS, SYMBOL_INDEX


def _dt(toy):
    ga = txtools.tx_load_bed(toy["bed"])
    genome = txtools.tx_load_genome(toy["fasta"])
    bam = txtools.tx_load_bam(toy["bam"], paired_end=True, load_seq=True, verbose=False)
    reads = txtools.tx_reads(bam, ga, min_reads=1, with_seq=True, verbose=False)
    dt = txtools.tx_makeDT_covNucFreq(reads, ga, genome=genome)
    return ga, genome, dt


def test_split_force_int():
    mat = np.zeros((3, len(IUPAC_CODE_2NUCS)), dtype=np.int64)
    # R=2 (even) -> 1 A + 1 G; Y=3 (odd) -> 1 C + 1 T + 1 N
    mat[0, SYMBOL_INDEX["R"]] = 2
    mat[1, SYMBOL_INDEX["Y"]] = 3
    out = simplify_force_int(mat)  # columns: A C G T N - .
    assert out[0, 0] == 1 and out[0, 2] == 1            # A, G
    assert out[1, 1] == 1 and out[1, 3] == 1 and out[1, 4] == 1  # C, T, N


def test_ratios_thresholds(toy):
    _, _, dt = _dt(toy)
    # default min_cov=50 with tiny coverage -> all NA
    dt2 = txtools.tx_add_startRatio(dt)
    assert dt2["startRatio"].isna().all()
    # min_cov=1 -> start_5p/cov
    dt3 = txtools.tx_add_startRatio(dt, min_cov=1)
    g1 = dt3[dt3.gene == "gene1"].reset_index(drop=True)
    assert abs(g1.loc[0, "startRatio"] - (2 / 2)) < 1e-9


def test_shifted_ratio_gene_boundaries(toy):
    _, _, dt = _dt(toy)
    dt = txtools.tx_add_startRatio1bpDS(dt, min_cov=1)
    # last position of each gene must be NA (nothing downstream within the gene)
    for g in ("gene1", "gene2"):
        sub = dt[dt.gene == g]
        assert np.isnan(sub["startRatio1bpDS"].iloc[-1])


def test_misinc(toy):
    _, _, dt = _dt(toy)
    dt = txtools.tx_add_misincRate(dt, min_nuc_reads=1, add_counts=True)
    row = dt[(dt.gene == "gene1") & (dt.txcoor == 5)].iloc[0]
    assert row["misincCount"] == 1
    assert row["nucTotal"] == 2
    assert abs(row["misincRate"] - 0.5) < 1e-9
    # a clean position has zero misincorporation
    clean = dt[(dt.gene == "gene1") & (dt.txcoor == 1)].iloc[0]
    assert clean["misincCount"] == 0


def test_gene_region(toy):
    ga, _, dt = _dt(toy)
    dt = txtools.tx_add_geneRegion(dt, ga)
    # thick spans the whole gene in the toy data -> all CDS
    assert (dt["geneRegion"] == "CDS").all()


def test_motif_presence(toy):
    ga, genome, dt = _dt(toy)
    dt = txtools.tx_add_motifPresence(dt, "CA", nuc_positions="all")
    col = "CA_motif_all"
    # every flagged run of length 2 should read 'CA' in refSeq for gene1
    g1 = dt[dt.gene == "gene1"].reset_index(drop=True)
    flagged = g1.index[g1[col]].tolist()
    seq = "".join(g1.refSeq.tolist())
    for i in flagged:
        # each flagged position is part of a CA occurrence
        assert seq[i] in "CA"


def test_rel_tx_pos(toy):
    _, _, dt = _dt(toy)
    dt = txtools.tx_add_relTxPos(dt)
    g1 = dt[dt.gene == "gene1"]
    assert abs(g1["relTxPos"].max() - 1.0) < 1e-9
    assert abs(g1[g1.txcoor == 30]["relTxPos"].iloc[0] - 0.5) < 1e-9
