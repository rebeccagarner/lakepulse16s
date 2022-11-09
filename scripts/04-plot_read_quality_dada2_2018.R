library(dada2)
packageVersion("dada2")


# Identify Cutadapt-trimmed samples and read files
(samples <- scan("samples", what = "character"))

(forward_reads <- paste0(samples, "_R1_trimmed.fastq.gz"))
(reverse_reads <- paste0(samples, "_R2_trimmed.fastq.gz"))


# Plot forward and reverse read quality
length_fwdreads <- length(forward_reads)
for (i in 1:round(length_fwdreads/25)) {
  pdf(paste0("lp2018watercolumn16s_fwdreads_quality_plots", (i*25-24), "-", (i*25), ".pdf"), width = 20, height = 20)
  print(plotQualityProfile(forward_reads[(i*25-24):(i*25)]))
  dev.off()
}

length_revreads <- length(reverse_reads)
for (i in 1:round(length_revreads/25)) {
  pdf(paste0("lp2018watercolumn16s_revreads_quality_plots", (i*25-24), "-", (i*25), ".pdf"), width = 20, height = 20)
  print(plotQualityProfile(reverse_reads[(i*25-24):(i*25)]))
  dev.off()
}
