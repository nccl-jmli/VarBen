

Format Introduction
==============================================================

### mutFile format
\#chrom  start       end         AF      type    alt
chr1	899778	    899778	    0.9	    snv	    T
chr1	1158637	    1158638	    0.9	    ins	    T
chr1	3397038	    3397039	    0.9	    ins	    T
chr1	3712508	    3712508	    0.9	    snv	    T
chr1	5927813	    5927813	    0.9	    snv	    T
chr1	6219495	    6219495	    0.9	    snv	    T
chr1	6504642	    6504642	    0.9	    snv	    T
chr1	6533124	    6533126	    0.9	    del	    .
chr1	7910956	    7910956	    0.9	    snv	    A
chr1	8399549	    8399549	    0.9	    snv	    A


### svFile format
#### del & inv format
\#chrom  start       end         type    AF
chrX	12994966	12996009	del     0.6
chrX	20172336	20176010	del     0.6
chrX	105121310	105134706	del	    0.6
chrX	108614726	108616334	del	    0.6
chrX    13703890    14134046    inv     0.6
chrX    19975999    20064786    inv     0.6
chrX    32391049    32794255    inv     0.6
chrX    40994338    41012689    inv     0.6

#### dup format
\#chrom  start       end         type    AF      dup_num
chr1    15808448    15814030    dup     0.6     3
chr1    16076907    16086182    dup     0.6     4
chr1    23665443    23711586    dup     0.6     3
chr1    28057278    28081157    dup     0.6     3

#### trans_chrom & trans_balance & trans_unbalance format
\#CHR1   CHR1_start  CHR1_end    type            AF      CHR2   CHR2_start   CHR2_end
chr10	7059511     7059511	    trans_chrom	    0.5	    chr19	17396810	17396810
chr19	17327977	17327977	trans_chrom     0.5	    chr3	186528041	186528041
chr3	107598967	107598967	trans_chrom	    0.5	    chr7	38371959	38371959
chr1	31561816	31561816	trans_chrom	    0.5	    chr6	41297838	41297838
chr2    29754284    29754947    trans_balance   0.5     chr2    42522695    42523089
chr10   43608984    43609308    trans_unbalance 0.5     chr6    117640981   117640982