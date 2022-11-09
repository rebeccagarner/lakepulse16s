library(dada2)
packageVersion("dada2")

# Identify Cutadapt-trimmed samples and read files
(samples <- scan("samples", what = "character"))

(forward_reads <- paste0(samples, "_R1_trimmed.fastq.gz"))
(reverse_reads <- paste0(samples, "_R2_trimmed.fastq.gz"))

# Prepare to filter reads
filtered_forward_reads <- paste0(samples, "_R1_filtered.fastq.gz")
filtered_reverse_reads <- paste0(samples, "_R2_filtered.fastq.gz")

# Define forward and reverse read lengths (based on quality plots)
forward_read_length <- 200
reverse_read_length <- 200

# Filter reads
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, maxEE = c(2, 2), rm.phix = TRUE, minLen = 160, truncLen = c(forward_read_length, reverse_read_length), compress = TRUE, multithread = TRUE)
