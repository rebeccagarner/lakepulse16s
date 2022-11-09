# LakePulse16S

Scripts associated with the bioinformatic analysis of

---

LakePulse 2017-2019 surface water 16S rRNA gene amplicon data (version *lp2017-2019watercolumn16sreseq2_636lakes*)

---

## Scripts

- **01-sample_list.R** Generate sample list(s).
- **02-concatenate_fastq.sh** Concatenate read files from the same samples.
- **03-cutadapt_2017/2018/2019.sh** Trim primer sequences from raw demultiplexed reads.
- **04-plot_read_quality_dada2_2017/2018/2019.R** Plot trimmed read quality.
- **05-filter_reads_dada2.R** Filter and trim reads from the 3' end.
- **06-dada2_2017/2018/2019.R** Infer amplicon sequence variants (ASVs) in dereplicated reads and merge forward and reverse paired-end reads.
- **07-combine_seqtabs_dada2.R** Combine ASV tables and remove bimeric ASVs.
- **08-count_reads_fastq.sh** Count number of reads in individual fastq files.
- **09-examine_nreads.R** Track read loss through the pipeline.
- **10a-taxass_wf.txt** Assign taxonomy (TaxAss workflow).
- **10b-reformat_dada2_seqtabs_modified.R** Format ASV table for TaxAss.
- **11-format_asvs.R** Format ASV and taxonomy tables.
- **12-align_asvs_mafft.txt** Construct a *de novo* multiple sequence alignment.
- **12-align_asvs_sina.txt** Align ASVs to SSU rRNA gene reference.
- **13-filter_asvs.R** Filter ASVs.
- **14-curate_samples.R** Filter ASVs based on taxonomy and filter samples.
