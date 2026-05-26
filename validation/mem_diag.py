#!/usr/bin/env python
"""Locate where Python txtools' peak memory is spent (per-phase RSS)."""
import sys, resource
import txtools
V = sys.argv[1]


def rss_gb():
    with open(f"/proc/{__import__('os').getpid()}/status") as fh:
        for ln in fh:
            if ln.startswith("VmRSS:"):
                return int(ln.split()[1]) / 1024 / 1024
    return float("nan")


def hwm_gb():
    with open(f"/proc/{__import__('os').getpid()}/status") as fh:
        for ln in fh:
            if ln.startswith("VmHWM:"):
                return int(ln.split()[1]) / 1024 / 1024
    return float("nan")


print(f"start                RSS={rss_gb():.2f} GB")
ga = txtools.tx_load_bed(f"{V}/chr19_refGene.bed12")
genome = txtools.tx_load_genome(f"{V}/chr19.fa")
bam = txtools.tx_load_bam(f"{V}/chr19.bam", paired_end=False, strand_mode=2,
                          load_seq=True, verbose=False)
print(f"after load handles   RSS={rss_gb():.2f} GB  (no BAM loaded into memory)")

reads = txtools.tx_reads(bam, ga, min_reads=50, with_seq=True, n_cores=8, verbose=False)
nfrag = reads.n_reads
seqbytes = sum(len(s) for g in reads for s in (g.seqs or []))
print(f"after tx_reads       RSS={rss_gb():.2f} GB  ({nfrag:,} fragments held, "
      f"~{seqbytes/1e9:.2f} GB of merged-seq characters)")

DT = txtools.tx_makeDT_covNucFreq(reads, ga, genome=genome)
print(f"after makeDT         RSS={rss_gb():.2f} GB  PEAK(VmHWM)={hwm_gb():.2f} GB")
mem = DT.memory_usage(deep=True)
print(f"final DataFrame:     {DT.shape[0]:,} rows x {DT.shape[1]} cols, "
      f"{mem.sum()/1e9:.2f} GB resident")
print("dtypes:", dict(DT.dtypes.astype(str)))
print(f"ru_maxrss peak       {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024/1024:.2f} GB")
