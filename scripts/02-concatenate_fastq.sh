#!/bin/bash
for sample in $(cat samples)
do
	echo "On sample: $sample"
	cat lp2019watercolumn16sreseq_raw/${sample}_R1.fastq.gz lp2019watercolumn16sreseq2_raw/${sample}_R1.fastq.gz > lp2019watercolumn16sreseq_concat_raw/${sample}_R1.fastq.gz
	cat lp2019watercolumn16sreseq_raw/${sample}_R2.fastq.gz lp2019watercolumn16sreseq2_raw/${sample}_R2.fastq.gz > lp2019watercolumn16sreseq_concat_raw/${sample}_R2.fastq.gz
done
