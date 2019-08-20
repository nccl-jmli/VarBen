#!/bin/bash


VARBEN_MUTEDITOR=../../bin/muteditor.py
SNV_INDEL_TEMPLATE=../input/chr21_snv_indel_template.tsv


python $VARBEN_MUTEDITOR -m $SNV_INDEL_TEMPLATE \
  -b ../input/IonTorrent_chr21.sorted.bam \
  -r ../reference/chr21.fa \
  -p 8 \
  --alignerIndex ../reference/chr21.fa \
  --seqer life \
  --aligner tmap \
  --haplosize 10 \
  --mindepth 500 \
  --minmutreads 10 \
  --snpfrac 0.2 \
  --libkey TCAG \
  --barcode TAAGGAGAAC \
  --floworder TACGTACGTCTGAGCATCGATCGATGTACAGC \
  -g \
  -o ../output/ion

