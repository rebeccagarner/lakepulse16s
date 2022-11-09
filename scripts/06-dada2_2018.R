library(dada2)
packageVersion("dada2")

# Identify Cutadapt-trimmed samples and read files
(samples <- scan("samples", what = "character"))

# Prepare to filter reads
filtered_forward_reads <- paste0(samples, "_R1_filtered.fastq.gz")
filtered_reverse_reads <- paste0(samples, "_R2_filtered.fastq.gz")

# Learn error rates
err_forward_reads <- learnErrors(filtered_forward_reads, nbases = 1e8, multithread = TRUE)
pdf("lp2018watercolumn16s_fwdreads_errors_plots.pdf")
plotErrors(err_forward_reads, nominalQ = TRUE)
dev.off()

err_reverse_reads <- learnErrors(filtered_reverse_reads, nbases = 1e8, multithread = TRUE)
pdf("lp2018watercolumn16s_revreads_errors_plots.pdf")
plotErrors(err_reverse_reads, nominalQ = TRUE)
dev.off()

# Dereplicate reads
derep_forward <- derepFastq(filtered_forward_reads, verbose = TRUE)
names(derep_forward) <- samples
derep_reverse <- derepFastq(filtered_reverse_reads, verbose = TRUE)
names(derep_reverse) <- samples

# Infer ASVs
dada_forward <- dada(derep_forward, err = err_forward_reads, pool = FALSE, multithread = TRUE)
dada_reverse <- dada(derep_reverse, err = err_reverse_reads, pool = FALSE, multithread = TRUE)

# Merge forward and reverse reads
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse, derep_reverse, trimOverhang = TRUE, minOverlap = 30)

# Create sequence table
seqtab <- makeSequenceTable(merged_amplicons)
dim(seqtab)
save(seqtab, file = "lp2018watercolumn16s_seqtab.rda")

# Track read loss through the pipeline
getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names = samples,
                          dada_fwd = sapply(dada_forward, getN), dada_rev = sapply(dada_reverse, getN),
                          merged = sapply(merged_amplicons, getN))
write.table(summary_tab, "lp2018watercolumn16s_dada2_readtracking.tsv", quote = FALSE, sep = "\t", col.names = NA)
