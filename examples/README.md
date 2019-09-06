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

Then go to `./output/ion` folder to check the results.

You will get three mutations in `edit.sorted.bam` file, one mutation is failed to spike in.

3. Spike SNV and Indel in Illumina BAM file

```
cd ./src
bash 03-spikein_snv_and_indel_Illumina.sh
```

4. Spike SV in WGS data

If you don't have any sample, please download NA12878 Low coverage WGS from [The International Genome Sample Resource](https://www.internationalgenome.org/data-portal/sample/NA12878). Then you should mapping the data to reference and create a index file for the BAM. After that, you could try to run the command below to edit SV in the BAM file.

```
cd ./src
bash 02-spikein_sv_Illumina_WGS.sh
```
