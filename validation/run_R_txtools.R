#!/usr/bin/env Rscript
# Reference run with the original txtools R package (1.0.6), matched params.
# txtools (pure R) comes from R_LIBS (a private lib holding only txtools);
# all compiled dependencies resolve from the consistent R 4.3.1 + Bioc 3.17 site.
suppressMessages(library(txtools))
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
V <- args[1]
nCores <- as.integer(args[2])
bam   <- file.path(V, "chr19.bam")
bed   <- file.path(V, "chr19_refGene.bed12")
fa    <- file.path(V, "chr19.fa")
outf  <- file.path(V, "r_chr19.tsv")

cat("txtools version:", as.character(packageVersion("txtools")), "\n")
t0 <- Sys.time()
geneAnnot <- tx_load_bed(bed)
genome    <- tx_load_genome(fa)
reads <- tx_load_bam(bam, pairedEnd = FALSE, loadSeq = TRUE,
                     strandMode = 2, loadSecondaryAligns = FALSE, verbose = TRUE)
cat("loaded", length(reads), "reads\n")

txReads <- tx_reads(reads, geneAnnot = geneAnnot, minReads = 50,
                    withSeq = TRUE, nCores = nCores, verbose = TRUE)
cat("genes with reads:", length(txReads), "\n")

DT <- tx_makeDT_covNucFreq(txReads, geneAnnot = geneAnnot, genome = genome,
                           nCores = nCores)
DT <- tx_add_misincRate(DT, minNucReads = 20, addCounts = TRUE)
DT <- tx_add_misincRateNucSpec(DT, refNuc = "C", misNuc = "T", minNucReads = 20)

fwrite(DT, outf, sep = "\t", na = "NA")
cat("rows:", nrow(DT), " cols:", ncol(DT), "\n")
cat("wrote", outf, "\n")
cat("R elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
