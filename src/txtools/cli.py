"""Command-line interface: the "simple to use executable".

Typical one-shot use::

    txtools makeDT --bam reads.bam --bed genes.bed12 --fasta genome.fa \\
        --table covNucFreq --paired --min-reads 50 --threads 8 \\
        --add misincRate,startRatio,geneRegion -o out.tsv

Run ``txtools makeDT --help`` for all options.
"""
from __future__ import annotations

import argparse
import sys
import time

from . import __version__
from . import api


# --------------------------------------------------------------------------- #
# --add metric dispatch
# --------------------------------------------------------------------------- #
def _apply_metric(df, spec: str, ctx: dict):
    """Apply one metric given a spec like ``name`` or ``name:arg1,arg2``."""
    if ":" in spec:
        name, rest = spec.split(":", 1)
        args = rest.split(",")
    else:
        name, args = spec, []
    name = name.strip()

    simple = {
        "startRatio": api.tx_add_startRatio,
        "endRatio": api.tx_add_endRatio,
        "startRatio1bpDS": api.tx_add_startRatio1bpDS,
        "startRatio1bpUS": api.tx_add_startRatio1bpUS,
        "endRatio1bpDS": api.tx_add_endRatio1bpDS,
        "endRatio1bpUS": api.tx_add_endRatio1bpUS,
        "nucTotal": api.tx_add_nucTotal,
        "misincCount": api.tx_add_misincCount,
        "misincRate": api.tx_add_misincRate,
        "relTxPos": api.tx_add_relTxPos,
        "pos": api.tx_add_pos,
    }
    if name in simple:
        return simple[name](df)
    if name == "geneRegion":
        return api.tx_add_geneRegion(df, ctx["gene_annot"])
    if name == "refSeqDT":
        return api.tx_add_refSeqDT(df, ctx["genome"], ctx["gene_annot"])
    if name == "misincRateNucSpec":
        if len(args) < 2:
            raise SystemExit("misincRateNucSpec needs ref,mis e.g. misincRateNucSpec:C,T")
        return api.tx_add_misincRateNucSpec(df, args[0], args[1])
    if name == "motifPresence":
        if not args:
            raise SystemExit("motifPresence needs a motif e.g. motifPresence:DRACH,center")
        motif = args[0]
        npos = args[1] if len(args) > 1 else "all"
        if npos not in ("all", "center"):
            npos = int(npos)
        return api.tx_add_motifPresence(df, motif, npos)
    if name == "rollingMean":
        if len(args) < 2:
            raise SystemExit("rollingMean needs col,win[,align] e.g. rollingMean:cov,50,center")
        align = args[2] if len(args) > 2 else "center"
        return api.tx_add_rollingMean(df, args[0], int(args[1]), align=align)
    raise SystemExit(f"Unknown metric in --add: {name!r}")


# --------------------------------------------------------------------------- #
# makeDT command
# --------------------------------------------------------------------------- #
def _cmd_makeDT(a: argparse.Namespace) -> int:
    t0 = time.time()
    paired = a.paired and not a.single
    if a.single:
        paired = False

    ga = api.tx_load_bed(a.bed)
    genome = api.tx_load_genome(a.fasta) if a.fasta else None
    need_seq = a.table in ("nucFreq", "covNucFreq")
    bam = api.tx_load_bam(
        a.bam, paired_end=paired, load_seq=need_seq,
        strand_mode=a.strand_mode, load_secondary=False, verbose=not a.quiet,
    )

    reads = api.tx_reads(
        bam, ga, min_reads=a.min_reads, with_seq=need_seq,
        ignore_strand=a.ignore_strand, n_cores=a.threads, verbose=not a.quiet,
    )

    if a.table == "coverage":
        df = api.tx_makeDT_coverage(reads, ga, genome=genome, full_dt=a.full)
    elif a.table == "nucFreq":
        df = api.tx_makeDT_nucFreq(reads, ga, genome=genome,
                                   simplify_IUPAC=a.simplify, full_dt=a.full)
    else:
        df = api.tx_makeDT_covNucFreq(reads, ga, genome=genome,
                                      simplify_IUPAC=a.simplify, full_dt=a.full)

    ctx = {"gene_annot": ga, "genome": genome}
    for spec in (a.add or []):
        spec = spec.strip()
        if spec:
            df = _apply_metric(df, spec, ctx)

    out = a.output
    if out in (None, "-"):
        df.to_csv(sys.stdout, sep="\t", index=False, na_rep="NA")
    else:
        df.to_csv(out, sep="\t", index=False, na_rep="NA")  # .gz inferred from name
        if not a.quiet:
            print(f"Wrote {len(df):,} rows x {df.shape[1]} cols to {out} "
                  f"in {time.time() - t0:.1f}s", file=sys.stderr)
    return 0


def _add_makeDT_parser(sub):
    p = sub.add_parser(
        "makeDT",
        help="BAM + BED12 + FASTA -> transcriptomic table (TSV).",
        description="Convert genomic alignments to a transcriptomic, "
                    "nucleotide-resolution table.",
    )
    req = p.add_argument_group("inputs")
    req.add_argument("--bam", required=True, help="Coordinate-sorted, indexed BAM.")
    req.add_argument("--bed", required=True, help="Gene annotation (BED-12 or BED-6).")
    req.add_argument("--fasta", help="Reference genome FASTA (adds refSeq column).")

    lib = p.add_argument_group("library")
    g = lib.add_mutually_exclusive_group()
    g.add_argument("--paired", action="store_true", default=True,
                   help="Paired-end reads (default).")
    g.add_argument("--single", action="store_true",
                   help="Single-end reads.")
    lib.add_argument("--strand-mode", type=int, default=1, choices=(1, 2),
                     dest="strand_mode",
                     help="1: strand from first mate (directional, default); "
                          "2: from last mate (dUTP/TruSeq).")
    lib.add_argument("--ignore-strand", action="store_true", dest="ignore_strand",
                     help="Assign reads to genes ignoring strand.")

    tab = p.add_argument_group("table")
    tab.add_argument("--table", default="covNucFreq",
                     choices=("coverage", "nucFreq", "covNucFreq"),
                     help="Which table to build (default: covNucFreq).")
    tab.add_argument("--min-reads", type=int, default=50, dest="min_reads",
                     help="Minimum fragments per gene (default: 50).")
    tab.add_argument("--simplify", default="splitForceInt",
                     choices=("splitForceInt", "splitHalf", "not"),
                     help="IUPAC ambiguity handling for nucleotide frequencies.")
    tab.add_argument("--full", action="store_true",
                     help="Include genes with no reads (zero rows).")
    tab.add_argument("--add", action="append", metavar="METRIC", default=None,
                     help="Add a derived metric column; repeatable. Examples: "
                          "--add misincRate --add startRatio --add geneRegion "
                          "--add misincRateNucSpec:C,T --add motifPresence:DRACH,center "
                          "--add rollingMean:cov,50,center")

    out = p.add_argument_group("output")
    out.add_argument("-o", "--output", help="Output TSV (.gz ok); '-' for stdout.")
    out.add_argument("-t", "--threads", type=int, default=1,
                     help="Worker processes (default: 1).")
    out.add_argument("-q", "--quiet", action="store_true", help="Suppress progress.")
    p.set_defaults(func=_cmd_makeDT)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="txtools",
        description="Transcriptomic, nucleotide-resolution analysis of RNA-seq "
                    "alignments (Python reimplementation of the txtools R package).",
    )
    p.add_argument("--version", action="version",
                   version=f"txtools {__version__}")
    sub = p.add_subparsers(dest="command", required=True)
    _add_makeDT_parser(sub)
    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
