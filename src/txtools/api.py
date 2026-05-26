"""Public API — the functions a user calls, with names kept consistent with
the original txtools R package.

    import txtools
    ga     = txtools.tx_load_bed("genes.bed12")
    genome = txtools.tx_load_genome("genome.fa")
    bam    = txtools.tx_load_bam("reads.bam", paired_end=True, load_seq=True)
    reads  = txtools.tx_reads(bam, ga, min_reads=50, with_seq=True, n_cores=8)
    DT     = txtools.tx_makeDT_covNucFreq(reads, ga, genome=genome)
    DT     = txtools.tx_add_misincRate(DT)
"""
from .io import tx_load_bam, tx_load_bed, tx_load_genome, Genome, GeneAnnotation, BamSource
from .core import tx_reads, TxReads, GeneReads
from .tables import (
    tx_makeDT_coverage,
    tx_makeDT_nucFreq,
    tx_makeDT_covNucFreq,
)
from .metrics import (
    tx_add_startRatio,
    tx_add_endRatio,
    tx_add_startRatio1bpDS,
    tx_add_startRatio1bpUS,
    tx_add_endRatio1bpDS,
    tx_add_endRatio1bpUS,
    tx_add_nucTotal,
    tx_add_misincCount,
    tx_add_misincRate,
    tx_add_misincRateNucSpec,
    tx_add_geneRegion,
    tx_add_motifPresence,
    tx_add_rollingMean,
    tx_add_relTxPos,
    tx_add_pos,
    tx_add_refSeqDT,
    tx_add_siteAnnotation,
)

__all__ = [
    "tx_load_bam", "tx_load_bed", "tx_load_genome",
    "Genome", "GeneAnnotation", "BamSource", "TxReads", "GeneReads",
    "tx_reads",
    "tx_makeDT_coverage", "tx_makeDT_nucFreq", "tx_makeDT_covNucFreq",
    "tx_add_startRatio", "tx_add_endRatio",
    "tx_add_startRatio1bpDS", "tx_add_startRatio1bpUS",
    "tx_add_endRatio1bpDS", "tx_add_endRatio1bpUS",
    "tx_add_nucTotal", "tx_add_misincCount", "tx_add_misincRate",
    "tx_add_misincRateNucSpec", "tx_add_geneRegion", "tx_add_motifPresence",
    "tx_add_rollingMean", "tx_add_relTxPos", "tx_add_pos",
    "tx_add_refSeqDT", "tx_add_siteAnnotation",
]
