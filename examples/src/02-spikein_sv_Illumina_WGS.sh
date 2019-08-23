#!/bin/bash


VARBEN_MUTEDITOR=../../bin/sveditor.py
SNV_INDEL_TEMPLATE=../input/wgs_sv_template.tsv



python $VARBEN_MUTEDITOR -m $SNV_INDEL_TEMPLATE \
	-b /home/yulijia/data/NA12878_low_coverage_WGS/wgs.sorted.bam \
	-r /home/admin/database/reference/hg19/ucsc.hg19.fasta \
	-p 8 \
	--alignerIndex /home/admin/database/reference/hg19/ucsc.hg19.fasta \
	--seqer illumina \
	--aligner bwa \
	--mindepth 2 \
	--minmutreads 2 \
	-l 101 \
	-o ../output/wgs
