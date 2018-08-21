import os
import argparse
from bameditor.deal_mut.checkMutInput import get_haplotypes
from bameditor.deal_mut.readsReplace import reads_replace
from bameditor.deal_mut.dealHaplotype import deal_haplotype_multi
from bameditor.deal_mut.readsModify import reads_modify
from bameditor.common.bamconvert import bamSort, bamMerge, bamIndex
from bameditor.common.recordlog import log


def main(run_args):
    invalid_log_file = os.path.join(run_args.outdir, 'invalid_mutation.txt')
    invalid_log = log(invalid_log_file)
    temp_out_dir = os.path.join(run_args.out_dir, "tempDir")
    os.system("mkdir -p %s" % temp_out_dir)

    # step1: deal with mutfile and get haplotypes
    haplotype_list = get_haplotypes(run_args.bamfile, run_args.reffasta, run_args.mutfile, int(run_args.haplosize),
                                    float(run_args.snpfrac), invalid_log)

    # step2: deal haplotypes and get total_chosen_reads, total_chosen_reads_muts
    total_chosen_reads, total_chosen_reads_muts = deal_haplotype_multi(run_args.bam_file, haplotype_list,
                                                                       temp_out_dir, run_args.reffasta,
                                                                       int(run_args.process), int(run_args.mindepth),
                                                                       int(run_args.minmutreads), int(run_args.minmapq),
                                                                       float(run_args.diffcover), run_args.single,
                                                                       run_args.is_multmapfilter, run_args.aligner,
                                                                       run_args.alignerIndex, invalid_log)

    # step3: modify the reads in total_chosen_reads itself
    reads_modify(total_chosen_reads, total_chosen_reads_muts, run_args.reffasta, int(run_args.prceoss))

    # step4: write edited reads to edited file and exclude reads to exclude file ,than remap edited file to reference
    edit_remap_bam_file, exclude_bam_file = reads_replace(run_args.bam_file, total_chosen_reads, run_args.seqer,
                                                          run_args.floworder, run_args.lib_key, run_args.barcode,
                                                          run_args.tag, temp_out_dir, run_args.aligner,
                                                          run_args.aligner_index, run_args.single)

    # step5: merge remap.edit.bam and exclude exclude.bam and sort
    bamIndex(exclude_bam_file)
    out_bam_file = os.path.join(temp_out_dir, "edit_exclude.bam")
    bamMerge([edit_remap_bam_file, exclude_bam_file], out_bam_file)

    out_sort_bam_file = os.path.join(run_args.outdir, "edit.sort.bam")
    out_sort_bam_file_prefix = os.path.join(run_args.outdir, "edit.sort")
    bamSort(out_bam_file, out_sort_bam_file_prefix)
    bamIndex(out_sort_bam_file)

    print "Edit Bam is completed! Result see %s and invalid mutation can't be spike in see %s." % (
        out_sort_bam_file, invalid_log_file)


def run():
    parser = argparse.ArgumentParser(description='Edit bamfile to spike in SNV, Indel, Complex, Substitution')
    parser.add_argument('-m', '--mutfile', dest='mutfile', required=True,
                        help='Target regions to try and spike in a mutation, format see README.txt')
    parser.add_argument('-b', '--bamfile', dest='bamfile', required=True, help='bam file')
    parser.add_argument('-r', '--reffasta', dest='reffasta', required=True, help='reference fasta file')
    parser.add_argument('-p', '--process', dest='process', required=False, default=1, help='process num(default = 1)')
    parser.add_argument('-g', '--single', action='store_true', default=False,
                        help="input BAM is simgle-ended (default is paired-end)")
    parser.add_argument('-o', '--outdir', dest='outdir', required=False, default="./",
                        help='out directory(default current path)')
    parser.add_argument('--seqer', dest='seqer', required=True, help='seqer (illumina, life, BGI)')
    parser.add_argument('--aligner', dest='aligner', required=True, default="bwa", help='choose a aligner(default bwa)')
    parser.add_argument('--alignerIndex', dest='alignerIndex', required=True, default=0.9, help='aligner genome index')

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
    parser.add_argument('--debug', action='store_true', default=False, help="debug mode(default False)")
    # parser.add_argument('--filetr_multimap', action='store_true', default=False, help="ignore reads with XA tag")
    # parser.add_argument('--filetr_mq', dest='filter_mq', default=30,
    #                     help="ignore reads with mq less than threshold (default 30)")

    run_args = parser.parse_args()
    assert os.path.exists(run_args.mutfile), "mutfile is not existed."
    assert os.path.exists(run_args.bamfile), "bamfile is not existed."
    assert os.path.exists(run_args.bamfile + ".bai"), "index of bam file is not existed."
    assert os.path.exists(run_args.reffasta), "reference fasta file is not existed."
    if run_args.debug:
        print type(run_args)
    main(run_args)

if __name__ == "__main__":
    run()
