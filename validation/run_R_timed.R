#!/usr/bin/env Rscript
# Per-phase timing of the R txtools pipeline (8 cores where supported).
suppressMessages(library(txtools)); suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE); V <- args[1]; nCores <- as.integer(args[2])
ph <- function(label, expr){ t <- Sys.time(); val <- force(expr)
  cat(sprintf("PHASE %-12s %7.1f s\n", label, as.numeric(difftime(Sys.time(), t, units="secs")))); val }

geneAnnot <- ph("load_bed",   tx_load_bed(file.path(V,"chr19_refGene.bed12")))
genome    <- ph("load_genome",tx_load_genome(file.path(V,"chr19.fa")))
reads     <- ph("load_bam",   tx_load_bam(file.path(V,"chr19.bam"), pairedEnd=FALSE,
                 loadSeq=TRUE, strandMode=2, loadSecondaryAligns=FALSE, verbose=FALSE))
txReads   <- ph("tx_reads",   tx_reads(reads, geneAnnot=geneAnnot, minReads=50,
                 withSeq=TRUE, nCores=nCores, verbose=FALSE))
DT        <- ph("makeDT",     tx_makeDT_covNucFreq(txReads, geneAnnot=geneAnnot,
                 genome=genome, nCores=nCores))
DT        <- ph("add_metrics",{ DT <- tx_add_misincRate(DT, minNucReads=20, addCounts=TRUE)
                 tx_add_misincRateNucSpec(DT, refNuc="C", misNuc="T", minNucReads=20) })
# Read THIS R process's own status (not a grep subshell's /proc/self).
status <- readLines(sprintf("/proc/%d/status", Sys.getpid()))
hwm <- grep("VmHWM", status, value = TRUE)
cat("PEAKMEM_SELF", hwm, "\n")  # peak resident set size of the R main process
cat("cores used:", nCores, "\n")
