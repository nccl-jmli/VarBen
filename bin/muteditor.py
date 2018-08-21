# Author: Shuangsang Fang
# Date: 2018/8/21
import os
import argparse
import sys
import time

scriptDir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..")
sys.path.append(scriptDir)

from bameditor.deal_mut.checkMutInput import get_haplotypes
from bameditor.deal_mut.readsReplace import reads_replace
from bameditor.deal_mut.dealHaplotype import deal_haplotype_multi
from bameditor.deal_mut.readsModify import reads_modify
from bameditor.common.bamconvert import bamSort, bamMerge, bamIndex, bamAddRG_picard, bamMerge_picard
from bameditor.common.recordlog import InvalidLog


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
                        help='Target regions to try and spike in a mutation, format see README.txt')
    parser.add_argument('-b', '--bamfile', dest='bamfile', required=True, help='bam file')
    parser.add_argument('-r', '--reffasta', dest='reffasta', required=True, help='reference fasta file')
    parser.add_argument('--seqer', dest='seqer', required=True, help='seqer (illumina, life, BGI)')
    parser.add_argument('--aligner', dest='aligner', required=True, default="bwa", help='choose a aligner(default bwa)')
    parser.add_argument('--alignerIndex', dest='alignerIndex', required=True, help='aligner genome index')
    parser.add_argument('-o', '--outdir', dest='outdir', required=True, help='out directory')

    parser.add_argument('-p', '--process', dest='process', required=False, default=1, help='process num(default = 1)')
    parser.add_argument('-g', '--single', action='store_true', default=False,
                        help="input BAM is single-ended (default is paired-end)")

    parser.add_argument('--haplosize', dest='haplosize', required=False, default=0,
                        help='haplotype size (default = 0)')
    parser.add_argument('--mindepth', dest='mindepth', required=False, default=30,
                        help='minimum depth(default = 1000)')
    parser.add_argument('--minmutreads', dest='minmutreads', required=False, default=5,
                        help='minimum depth(default = 5)')
    parser.add_argument('--snpfrac', dest='snpfrac', required=False, default=1, help='snp fraction(default = 1)')
    parser.add_argument('--minmapq', dest='minmapq', required=False, default=20,
                        help='read mapping quality less than minmapq will not be considered to edit')
    parser.add_argument('--multmapfilter', action='store_true', default=False,
                        help="multiple mapped reads will not be considered to edit. (default is True)")
    parser.add_argument('--floworder', dest='floworder', required=False, help='flower order of life sequence')
    parser.add_argument('--libkey', dest='libkey', required=False, help='libkey of life sequence')
    parser.add_argument('--barcode', dest='barcode', required=False, help='barcode of life sequence')
    parser.add_argument('--diffcover', dest='diffcover', required=False, default=0.6,
                        help='coverage difference(default 0.9)')
    parser.add_argument('--tag', action='store_true', default=False, help="tag mutated reads (default False)")
    # parser.add_argument('--debug', action='store_true', default=False, help="debug mode(default False)")

    run_args = parser.parse_args()
    assert os.path.exists(run_args.mutfile), "mutfile is not existed."
    assert os.path.exists(run_args.bamfile), "bamfile is not existed."
    assert os.path.exists(run_args.bamfile + ".bai"), "index of bam file is not existed."
    assert os.path.exists(run_args.reffasta), "reference fasta file is not existed."

    main(run_args)


if __name__ == "__main__":
    run()
