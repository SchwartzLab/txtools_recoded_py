#!/usr/bin/env python
"""Per-phase timing of the Python txtools pipeline (matched to the R phases)."""
import sys, time
import txtools
V, ncores = sys.argv[1], int(sys.argv[2])


def ph(label, fn):
    t = time.time(); val = fn(); print(f"PHASE {label:<12} {time.time()-t:7.1f} s", flush=True); return val


ga = ph("load_bed", lambda: txtools.tx_load_bed(f"{V}/chr19_refGene.bed12"))
genome = ph("load_genome", lambda: txtools.tx_load_genome(f"{V}/chr19.fa"))
bam = ph("load_bam", lambda: txtools.tx_load_bam(f"{V}/chr19.bam", paired_end=False,
         load_seq=True, strand_mode=2, verbose=False))
reads = ph("tx_reads", lambda: txtools.tx_reads(bam, ga, min_reads=50, with_seq=True,
           n_cores=ncores, verbose=False))
DT = ph("makeDT", lambda: txtools.tx_makeDT_covNucFreq(reads, ga, genome=genome))


def _add(df):
    df = txtools.tx_add_misincRate(df, min_nuc_reads=20, add_counts=True)
    return txtools.tx_add_misincRateNucSpec(df, "C", "T", min_nuc_reads=20)


DT = ph("add_metrics", lambda: _add(DT))
import resource
self_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024
child_gb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024
print(f"PEAKMEM_SELF {self_gb:.3f} GB")
print(f"PEAKMEM_CHILD_MAX {child_gb:.3f} GB")
print("cores used:", ncores)

