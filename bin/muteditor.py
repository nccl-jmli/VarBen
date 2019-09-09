# Author: Shuangsang Fang
# Date: 2018/8/21
import os
import argparse
import sys
import time

scriptDir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..")
sys.path.append(scriptDir)

from varben.deal_mut.checkMutInput import get_haplotypes
from varben.deal_mut.readsReplace import reads_replace
from varben.deal_mut.dealHaplotype import deal_haplotype_multi
from varben.deal_mut.readsModify import reads_modify
from varben.common.bamconvert import bamSort, bamMerge, bamIndex, bamAddRG_picard, bamMerge_picard
from varben.common.recordlog import InvalidLog


def main(run_args):
    start_time = time.asctime(time.localtime(time.time()))
    # print start_time
    temp_out_dir = os.path.join(run_args.outdir, "tempDir")
    os.system("mkdir -p %s" % temp_out_dir)
    invalid_log_file = os.path.join(run_args.outdir, 'invalid_mutation.txt')
    invalid_log = InvalidLog(invalid_log_file)

    # step1: deal with mutfile and get haplotypes
    print "step1: deal with mutfile and get haplotypes"
    haplotype_list = get_haplotypes(run_args.bamfile, run_args.reffasta, run_args.mutfile, int(run_args.haplosize),
                                    float(run_args.snpfrac), invalid_log)

    # step2: deal haplotypes and get total_chosen_reads, total_chosen_reads_muts
    print "step2: deal haplotypes and get total_chosen_reads, total_chosen_reads_muts"
    success_list_file = os.path.join(run_args.outdir, 'success_list.txt')
    total_chosen_reads, total_chosen_reads_muts = deal_haplotype_multi(run_args.bamfile, haplotype_list,
                                                                       temp_out_dir, run_args.reffasta,
                                                                       int(run_args.process),
                                                                       int(run_args.mindepth),
                                                                       int(run_args.minmutreads),
                                                                       int(run_args.minmapq),
                                                                       float(run_args.diffcover),
                                                                       run_args.single,
                                                                       run_args.multmapfilter, run_args.aligner,
                                                                       run_args.alignerIndex, invalid_log,
                                                                       success_list_file)
    invalid_log.close()
    if len(total_chosen_reads) == 0:
        print "Warning: No reads to deal with of all these sv, checkout your sv file"
        return

    # step3: modify the reads in total_chosen_reads itself
    print "step3: modify the reads in total_chosen_reads itself"
    reads_modify(total_chosen_reads, total_chosen_reads_muts, run_args.reffasta, int(run_args.process))

    # step4: write edited reads to edited file and exclude reads to exclude file ,than remap edited file to reference
    print "step4: write edited reads to edited file and exclude reads to exclude file ,than remap edited file to reference"
    edit_remap_bam_file, exclude_bam_file = reads_replace(run_args.bamfile, total_chosen_reads, run_args.seqer,
                                                          run_args.floworder, run_args.libkey, run_args.barcode,
                                                          run_args.tag, temp_out_dir, run_args.aligner,
                                                          run_args.alignerIndex, run_args.single)

    # step5: merge remap.edit.bam and exclude exclude.bam and sort
    print "step5: merge remap.edit.bam and exclude exclude.bam and sort"
    # edit_remap_bam_file, exclude_bam_file = os.path.join(temp_out_dir, "edit.remap.sort.bam"), os.path.join(
    #     temp_out_dir, "exclude.bam")
    out_bam_file = os.path.join(run_args.outdir, "edit.sorted.bam")
    bamMerge([edit_remap_bam_file, exclude_bam_file], out_bam_file)
    bamIndex(out_bam_file)
    end_time = time.asctime(time.localtime(time.time()))
    # speed_time = end_time - start_time
    print "Edit Bam is completed! Result see %s and valid mutation see %s. Invalid mutation can't be spike in see %s." % (
        out_bam_file, success_list_file, invalid_log_file)


def run():
    parser = argparse.ArgumentParser(description='Edit bamfile to spike in SNV, Indel, Complex, Substitution')
    parser.add_argument('-m', '--mutfile', dest='mutfile', required=True,
                        help='Target regions to try and spike in a point mutation, format see README.txt')
    parser.add_argument('-b', '--bamfile', dest='bamfile', required=True,
                        help='The input BAM files should be aligned to reference genome, sorted in coordinate order and indexed with samtools index. By default, the software considers the BAM file consists by entirely paired-end reads, if a user needs to spike mutation in a BAM file which consists by single-end reads, they need using the -single option.')
    parser.add_argument('-r', '--reffasta', dest='reffasta', required=True, help='The reference genome in FASTA format and indexed with samtools (samtools faidx). he target BAM file should be generated by the same reference file used in this option, especially the chromosome names and lengths in the reference FASTA must be the same as in the BAM header.')
    parser.add_argument('-o', '--outdir', dest='outdir', required=True, help='A output directory name for edited bam file and other information.')
    parser.add_argument('--alignerIndex', dest='alignerIndex', required=True, help='The indexed reference genome in the FASTA format of aligner. For example, if the aligner is BWA, then BWA index should be provided. This FASTA file is called by the external aligner.')
    parser.add_argument('-p', '--process', dest='process', required=False, default=1,
                        help='Parallel mode: process number (default = 1).')

    parser.add_argument('--seqer', dest='seqer', required=False, default="illumina", help='Define the sequence platform for reads generating: illumina, life, BGI. Default is illumina, life stand for the Ion Torrent platform.')
    parser.add_argument('-g', '--single', action='store_true', default=False,
                        help="To declare that the input bam is single-ended (default is False)")
    parser.add_argument('--aligner', dest='aligner', required=False, default="bwa",
                        help='Choose an aligner from bwa, novoalign and tmap (default bwa)')

    parser.add_argument('--haplosize', dest='haplosize', required=False, default=0,
                        help='The size of haplotype block to consider when adding more than 1 proximal mutation. (default = 0) For example, if two SNVs are spiked in 5bp apart and -haplosize is 5 or greater, the two SNVs will be on the same haplotype (i.e. share the same reads for reads covering both positions).')
    parser.add_argument('--mindepth', dest='mindepth', required=False, default=30,
                        help='The minimum depth of reads position which could be edited to simulate mutation. (default = 30)')
    parser.add_argument('--minmutreads', dest='minmutreads', required=False, default=5,
                        help='The minimum number of reads to be edited in one position (default = 5). VarBen will calculate the number of mutated reads by the allele frequency and the total number of reads in the position. If the mutation reads number is less than 5 and --minmutreads is 5 or greater, VarBen will drop this site automatically.')
    parser.add_argument('--snpfrac', dest='snpfrac', required=False, default=1, help='To avoid spike any mutatoin on top of existing heterozygous alleles, the heterozygous allele fraction set to 0.1 (default = 1)')
    parser.add_argument('--minmapq', dest='minmapq', required=False, default=20,
                        help='Reads mapping quality less than MINMAPQ will not be considered to edit (default 20).')
    parser.add_argument('--multmapfilter', action='store_true', default=False,
                        help="Any multi-mapped reads will not be considered to edit (default is True).")
    parser.add_argument('--diffcover', dest='diffcover', required=False, default=0.6,
                        help='The coverage difference allowed between the input BAM and output BAM (default 0.9).')
    parser.add_argument('--floworder', dest='floworder', required=False, help='If the sequence platform is Ion torrent (--seqer is life), sequencing flower order should be provided.')
    parser.add_argument('--libkey', dest='libkey', required=False, help='If the sequence platform is Ion torrent (--seqer is life), the library key sequence should be provided.')
    parser.add_argument('--barcode', dest='barcode', required=False, help='If the sequence platform is Ion torrent (--seqer is life), the library barcode sequence should be provided.')
    parser.add_argument('--tag', action='store_true', default=False, help="Add tag to edited reads (default False).")
    # parser.add_argument('--debug', action='store_true', default=False, help="debug mode(default False)")

    run_args = parser.parse_args()
    assert os.path.exists(run_args.mutfile), "mutfile is not existed."
    assert os.path.exists(run_args.bamfile), "bamfile is not existed."
    assert os.path.exists(run_args.bamfile + ".bai"), "index of bam file is not existed."
    assert os.path.exists(run_args.reffasta), "reference fasta file is not existed."

    main(run_args)


if __name__ == "__main__":
    run()
