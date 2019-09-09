# Author: Shuangsang Fang
# Date: 2018/8/21
import os
import argparse
import time
import sys

scriptDir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..")
sys.path.append(scriptDir)
from varben.common.methods import get_insertSize_range
from varben.deal_sv.dealSVType import deal_sv
from varben.common.bamconvert import bamIndex, bamMerge
from varben.deal_sv.checkSVInput import check_sv_file
from varben.common.recordlog import InvalidLog, RunLog
from varben.deal_sv.getReadsByRegion import get_reads_by_region
from varben.deal_sv.writeBamByChr import write_sub_bam
from varben.deal_sv.mergeEditBam import merge_edit_bam


def main(run_args):
    start_time = time.asctime(time.localtime(time.time()))
    # print start_time
    if not os.path.exists(run_args.outdir):
        os.mkdir(run_args.outdir)
    invalid_log_file = os.path.join(run_args.outdir, 'invalid_mutation.txt')
    invalid_log = InvalidLog(invalid_log_file)

    run_log_file = os.path.join(run_args.outdir, 'run.log')
    run_log = RunLog(run_log_file)

    temp_out_dir = os.path.join(run_args.outdir, "tempDir")
    if not os.path.exists(temp_out_dir):
        os.mkdir(temp_out_dir)

    # step0: prepare sv list
    sv_list = check_sv_file(run_args.svfile, run_args.reffasta, invalid_log)

    if not sv_list:
        exit("no sv list to deal with")

    # step1: get insert size of paired reads
    print "step1: get insert size of paired reads"
    insert_size = get_insertSize_range(run_args.bamfile, run_args.readlength, run_args.single)

    # step2: deal with sv and get total edited reads
    print "step2: deal with sv and get total edited reads"
    success_file = os.path.join(run_args.outdir, 'success_list.txt')
    total_modify_reads_file, total_delete_reads_file, total_add_reads_file, total_modify_list, total_delete_list, total_add_list = deal_sv(
        run_args.bamfile, run_args.reffasta, sv_list,
        run_args.single,
        int(run_args.minmapq),
        run_args.multmapfilter,
        int(run_args.mindepth),
        int(run_args.minmutreads),
        int(run_args.readlength),
        temp_out_dir,
        insert_size,
        invalid_log,
        run_log, success_file)
    invalid_log.close()

    # step3: get reads by region bed and write bam file
    print "step3: get reads by region bed and write bam file"
    chrom_list, used_bam_file_tmp, exclude_bam_file_tmp = get_reads_by_region(run_args.bamfile, sv_list,
                                                                              temp_out_dir)

    # write reads which may probably used to used.bam and reads should not be used to exclude.bam
    used_bam_file, exclude_bam_file = write_sub_bam(chrom_list, used_bam_file_tmp, exclude_bam_file_tmp,
                                                    temp_out_dir, total_modify_reads_file,
                                                    total_delete_reads_file, total_add_reads_file,
                                                    int(run_args.process))

    # step4: merge edited reads and remap to new bam, consider about the tag, RG, life reads
    print "step4: merge edited reads and remap to new bam, consider about the tag, RG, life reads"
    edit_remap_bam_file = merge_edit_bam(run_args.bamfile, temp_out_dir, run_args.single, total_modify_reads_file,
                                         total_add_reads_file, used_bam_file, total_modify_list, total_add_list,
                                         run_args.seqer, run_args.aligner, run_args.alignerIndex,
                                         run_args.floworder,
                                         run_args.libkey, run_args.barcode, run_args.tag)

    # step5: remapped edit reads and merge
    print "step5: remapped edit reads and merge"
    out_bam_file = os.path.join(run_args.outdir, "edit.sorted.bam")
    bamMerge([edit_remap_bam_file, exclude_bam_file], out_bam_file)
    bamIndex(out_bam_file)

    end_time = time.asctime(time.localtime(time.time()))
    # print end_time
    # speed_time = end_time - start_time
    print "Edit Bam is completed! Result see %s and valid mutation see %s. Invalid mutation can't be spike in see %s." % (
        out_bam_file, success_file, invalid_log_file)


def run():
    parser = argparse.ArgumentParser(description='Edit bam file to spike in SV')
    parser.add_argument('-m', '--svfile', dest='svfile', required=True,
                        help='Target regions to try and spike in a SV, format see README.txt')
    parser.add_argument('-b', '--bamfile', dest='bamfile', required=True, help='The input BAM files should be aligned to reference genome, sorted in coordinate order and indexed with samtools index. By default, the software considers the BAM file consists by entirely paired-end reads, if a user needs to spike mutation in a BAM file which consists by single-end reads, they need using the -single option.')
    parser.add_argument('-r', '--reffasta', dest='reffasta', required=True, help='The reference genome in FASTA format and indexed with samtools (samtools faidx). he target BAM file should be generated by the same reference file used in this option, especially the chromosome names and lengths in the reference FASTA must be the same as in the BAM header.')
    parser.add_argument('-l', '--readlength', dest='readlength', required=True, help='read length of bam reads')
    parser.add_argument('-o', '--outdir', dest='outdir', required=True, help='A output directory name for edited bam file and other information.')
    parser.add_argument('--alignerIndex', dest='alignerIndex', required=True, help='The indexed reference genome in the FASTA format of aligner. For example, if the aligner is BWA, then BWA index should be provided. This FASTA file is called by the external aligner.')
    parser.add_argument('-p', '--process', dest='process', required=False, default=1,
                        help='Parallel mode: process number (default = 1).')
    parser.add_argument('--seqer', dest='seqer', default='illumina', help='Define the sequence platform for reads generating: illumina, life, BGI. Default is illumina, life stand for the Ion Torrent platform.')
    parser.add_argument('-g', '--single', action='store_true', default=False,
                        help="To declare that the input bam is single-ended (default is False)")
    parser.add_argument('--aligner', dest='aligner', default="bwa",
                        help='Choose an aligner from bwa, novoalign and tmap (default bwa)')

    parser.add_argument('--mindepth', dest='mindepth', required=False, default=30,
                        help='The minimum depth of reads position which could be edited to simulate mutation. (default = 1000)')
    parser.add_argument('--minmutreads', dest='minmutreads', required=False, default=5,
                        help='The minimum number of reads to be edited in one position (default = 5). VarBen will calculate the number of mutated reads by the allele frequency and the total number of reads in the position. If the mutation reads number is less than 5 and --minmutreads is 5 or greater, VarBen will drop this site automatically.')
    parser.add_argument('--minmapq', dest='minmapq', required=False, default=20,
                        help='Reads mapping quality less than MINMAPQ will not be considered to edit (default 20).')
    parser.add_argument('--multmapfilter', action='store_true', default=False,
                        help="Any multi-mapped reads will not be considered to edit (default is True).")

    parser.add_argument('--floworder', dest='floworder', required=False, help='If the sequence platform is Ion torrent (--seqer is life), sequencing flower order should be provided.')
    parser.add_argument('--libkey', dest='libkey', required=False, help='If the sequence platform is Ion torrent (--seqer is life), the library key sequence should be provided.')
    parser.add_argument('--barcode', dest='barcode', required=False, help='If the sequence platform is Ion torrent (--seqer is life), the library barcode sequence should be provided.')
    parser.add_argument('--tag', action='store_true', default=False, help="Add tag to edited reads (default False).")
    run_args = parser.parse_args()
    assert os.path.exists(run_args.svfile), "svfile is not existed."
    assert os.path.exists(run_args.bamfile), "bamfile is not existed."
    assert os.path.exists(run_args.bamfile + ".bai") or os.path.exists(run_args.bamfile[:run_args.bamfile.rindex(".")]+".bai"), "index of bam file is not existed."
    assert os.path.exists(run_args.reffasta), "reference fasta file is not existed."

    main(run_args)


if __name__ == "__main__":
    run()
