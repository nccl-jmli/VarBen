Examples
=============

1. please first to unzip all reference file in ./reference folder

```
xz -d ./reference/chr21*
```

2. Spike SNV and Indel in Ion Torrent BAM file (target sequencing)

```
cd ./src
bash 01-spikein_snv_and_indel_IonTorrent.sh
```

Then go to ./output/ion folder to check the results.

You will get three mutations in `edit.sorted.bam` file, one mutation is failed to spike in (due to the panel didn't cover this region) 

