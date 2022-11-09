#!/bin/bash
# Primer 515F: GTGCCAGCMGCCGCGGTAA
# Primer 806R: GGACTACHVGGGTWTCTAAT
cutadapt --version  # 3.7
for sample in $(cat samples)
do
	echo "On sample: $sample"
	cutadapt -a GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC -A GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC -m 160 -M 231 --discard-untrimmed -o lp2019watercolumn16sreseq_concat_trimmed/${sample}_R1_trimmed.fastq.gz -p lp2019watercolumn16sreseq_concat_trimmed/${sample}_R2_trimmed.fastq.gz ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz >> lp2019watercolumn16sreseq_concat_trimmed/lp2019watercolumn16sreseq_concat_cutadapt_stats.txt 2>&1
done

# Summarize the fraction of reads and base pairs retained in each sample
paste <(grep "passing" lp2019watercolumn16sreseq_concat_trimmed/lp2019watercolumn16sreseq_concat_cutadapt_stats.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" lp2019watercolumn16sreseq_concat_trimmed/lp2019watercolumn16sreseq_concat_cutadapt_stats.txt | cut -f3 -d "(" | tr -d ")") > lp2019watercolumn16sreseq_concat_trimmed/lp2019watercolumn16sreseq_concat_cutadapt_summary.txt

# In R:
# library(tidyverse)
# 
# samples <- read_tsv("samples", col_names = "sample_id")
# 
# summary <- read_tsv("lp2019watercolumn16sreseq_concat_cutadapt_summary.txt", col_names = c("pct_passing", "pct_filtered")) %>%
#   mutate(pct_passing = as.numeric(str_remove(pct_passing, "%")),
#          pct_filtered = as.numeric(str_remove(pct_filtered, "%")))
# 
# sample_summary <- samples %>%
#   bind_cols(summary)
# 
# sample_summary %>%
#   write_tsv("lp2019watercolumn16sreseq_concat_cutadapt_summary.txt", col_names = TRUE)