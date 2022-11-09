#!/bin/bash
for filename in *.fastq.gz
do
	prefix=`basename $filename .fastq.gz`
	zcat $filename | echo $((`wc -l`/4)) > lp2017-2019watercolumn16s_nreads/${prefix}_nreads.txt
done
