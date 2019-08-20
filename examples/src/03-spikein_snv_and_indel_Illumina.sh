#!/bin/bash

VARBEN_MUTEDITOR=../../bin/muteditor.py
SNV_INDEL_TEMPLATE=../input/EGFR.tsv

python $VARBEN_MUTEDITOR -m $SNV_INDEL_TEMPLATE \
	-b ../input/EGFR_sort.bam \
	-r ../reference/chr7.fa \
	-p 2 \
	--alignerIndex ../reference/chr7.fa \
	--seqer illumina \
	--aligner bwa \
	--haplosize 10 \
	--mindepth 100 \
	--minmutreads 10 \
	--snpfrac 0.1 \
	-o ../output/illumina
