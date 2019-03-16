VarBen
=========================
A tool for adding variant simulation to .bam files, including single-nucleotide variants, short insertions and deletions, and large structural variants.

As for SV editing, we can use VarBen to add variations like deletion, inversion, duplication, translocation (include balance translocation, unbalance translocation, chromosome translocation) in .bam files.


Installation
-------------------------
Dependencies: 
+ pysam (python package >=0.9.4)

```
pip install pysam
```

+ samtools (http://samtools.sourceforge.net/)

```
git clone https://github.com/samtools/htslib.git
make -C htslib

git clone https://github.com/samtools/samtools.git
make -C samtools
cp samtools/samtools $HOME/bin

git clone https://github.com/samtools/bcftools.git
make -C bcftools
cp bcftools/bcftools $HOME/bin
```

+ bwa (http://bio-bwa.sourceforge.net/)

```
git clone https://github.com/lh3/bwa.git
make -C bwa
cp bwa/bwa $HOME/bin
```

Quick Start
---------------------
+ Spike point mutations (SNV & InDel) into bam file. (Illumina platform)

```bash
python /opt/VarBen/bin/muteditor.py -m ./mutFile.tsv \
-b ./Illumina_normal.bam \
-r ./reference/hg19/ucsc.hg19.fasta \
-p 4 \
--aligner bwa \
--alignerIndex ./reference/hg19/ucsc.hg19.fasta \
--seqer illumina \
--haplosize 10 \
--mindepth 500 \
--minmutreads 5 \
--snpfrac 0.1 \
-o ./Illumina_mut_out/

```
+ Spike point mutations (SNV & InDel) into bam file. (Ion torrent platform)

```bash
python /opt/VarBen/bin/muteditor.py -m ./mutFile.tsv \
-b ./IonXpress_001_realigned.bam \
-r ./referenceLibrary/tmap-f3/hg19/hg19.fasta \
-p 4 \
--aligner tmp \
--alignerIndex ./referenceLibrary/tmap-f3/hg19/hg19.fasta \
--seqer life \
--haplosize 10 \
--mindepth 500 \
--minmutreads 5 \
--snpfrac 0.1 \
--libkey TCAG \
--barcode CTAAGGTAACGAT \
--floworder TACGTACGTCTGAGCATCGATCGATGTACAGC \
-g \
-o ./Ion_mut_out/

```

+ Spike SVs into bam file.  (Illumina platform)

```bash
python /opt/VarBen/bin/sveditor.py -m ./mutFile.tsv \
-b ./Illumina_normal.bam \
-r ./reference/hg19/ucsc.hg19.fasta \
--aligner bwa \
--alignerIndex ./reference/hg19/ucsc.hg19.fasta \
--seqer illumina \
--mindepth 100 \
--minmutreads 10 \
-p 12 \
-l 100 \
-o ./sv_out/

```


Format Introduction
-------------------------
### mutFile format

```
#chrom  start end AF  type  alt

chr1  899778  899778  0.9 snv T

chr1  3712508 3712508 0.9 snv T

chr1  1158637 1158638 0.9 ins TAG

chr1	3397038 3397039 0.9 ins AGGTAG

chr1	6533124 6533126 0.9 del .

chr1	7910946 7910956 0.9 del .
chr7	55242467	55242481	0.3	Sub	TTC	
```

### svFile format
#### del & inv format

```
#chrom  start end type  AF

chrX  12994966  12996009  del 0.6

chrX  20172336  20176010  del 0.6

chrX  105121310	105134706	del 0.6

chrX  108614726	108616334	del 0.6

chrX  13703890  14134046  inv 0.6

chrX  19975999  20064786  inv 0.6

chrX  32391049  32794255  inv 0.6

chrX  40994338  41012689  inv 0.6
```


#### dup format

```
#chrom  start end type  AF  dup_num

chr1  15808448  15814030  dup 0.6 3

chr1  16076907  16086182  dup 0.6 4

chr1  23665443  23711586  dup 0.6 3

chr1  28057278  28081157  dup 0.6 3
```

#### trans_chrom & trans_balance & trans_unbalance format

```
#CHR1 CHR1_start  CHR1_end  type  AF  CHR2  CHR2_start  CHR2_end

chr10 7059511 7059511 trans_chrom 0.5 chr19 17396810  17396810

chr19 17327977  17327977  trans_chrom 0.5 chr3  186528041 186528041

chr3  107598967 107598967 trans_chrom 0.5 chr7  38371959  38371959

chr1  31561816  31561816  trans_chrom 0.5 chr6  41297838  41297838

chr2  29754284  29754947  trans_balance 0.5 chr2  42522695  42523089

chr10 43608984  43609308  trans_unbalance 0.5 chr6  117640981 117640982
```

#### CNV format

```
#chrom  start end type  AF  cnv_type

chrX 66764255 66950650 cnv 2.5 gain 

chr20 52186265 52200826 cnv 2 lose

```

