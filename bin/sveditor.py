import os
import argparse
import time
import sys

scriptDir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..")
sys.path.append(scriptDir)
from bameditor.common.methods import get_insertSize_range
from bameditor.deal_sv.dealSVType import deal_sv
from bameditor.common.bamconvert import bamIndex, bamMerge_picard, bamAddRG_picard
from bameditor.deal_sv.checkSVInput import check_sv_file
from bameditor.common.recordlog import InvalidLog, RunLog
from bameditor.deal_sv.getReadsByRegion import get_reads_by_region
from bameditor.deal_sv.writeBamByChr import write_sub_bam
from bameditor.deal_sv.mergeEditBam import merge_edit_bam


def main(run_args):
    start_time = time.asctime(time.localtime(time.time()))
    print start_time
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
    if run_args.debug:
        print len(sv_list)
    if not sv_list:
        exit("no sv list to deal with")
    
    # step1: get insert size of paired reads
    print "step1: get insert size of paired reads"
    insert_size = get_insertSize_range(run_args.bamfile, run_args.readlength, run_args.single, run_args.debug)
    if run_args.debug:
        print insert_size

    # step2: deal with sv
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
    if run_args.debug:
        print "total modify/add/delete reads file", total_modify_reads_file, total_delete_reads_file, total_add_reads_file

    # step3: get reads by region bed and write bam file
    print "step3: get reads by region bed and write bam file"
    chrom_list, used_bam_file_tmp, exclude_bam_file_tmp = get_reads_by_region(run_args.bamfile, sv_list,
                                                                              temp_out_dir)
    
    # write reads which may probably used to used.bam and reads should not be used to exclude.bam
    used_bam_file, exclude_bam_file = write_sub_bam(chrom_list, used_bam_file_tmp, exclude_bam_file_tmp,
                                                    temp_out_dir, total_modify_reads_file,
                                                    total_delete_reads_file, total_add_reads_file,
                                                    int(run_args.process))

    if run_args.debug:
        print "used & exclude bam:", used_bam_file, exclude_bam_file
    
    # step4: merge edited reads and remap to new bam, consider about the tag, RG, life reads
    print "step4: merge edited reads and remap to new bam, consider about the tag, RG, life reads"
    edit_remap_bam_file = merge_edit_bam(run_args.bamfile, temp_out_dir, run_args.single, total_modify_reads_file,
                                         total_add_reads_file, used_bam_file, total_modify_list, total_add_list,
                                         run_args.seqer, run_args.aligner, run_args.alignerIndex,
                                         run_args.floworder,
                                         run_args.libkey, run_args.barcode, run_args.tag)

    if run_args.debug:
        print "edit remap bam:", edit_remap_bam_file

    # step5: remapped edit reads and merge
    print "step5: remapped edit reads and merge"
    edit_remap_addRG_bam_file = os.path.join(temp_out_dir, "edit.remap.addRG.bam")
    bamAddRG_picard(edit_remap_bam_file, edit_remap_addRG_bam_file, run_args.picard_path, run_args.seqer)
    out_bam_file = os.path.join(run_args.outdir, "edit.sorted.bam")
    bamMerge_picard([edit_remap_addRG_bam_file, exclude_bam_file], out_bam_file, run_args.picard_path,
                    sort_order='coordinate')
    bamIndex(out_bam_file)

    end_time = time.asctime(time.localtime(time.time()))
    print end_time


def run():
    parser = argparse.ArgumentParser(description='Edit bamfile to spike in SNV, Indel, Complex, Substitution')
    parser.add_argument('-m', '--svfile', dest='svfile', required=True,
                        help='Target regions to try and spike in a sv, format see README.txt')
    parser.add_argument('-b', '--bamfile', dest='bamfile', required=True, help='bam file')
    parser.add_argument('-r', '--reffasta', dest='reffasta', required=True, help='reference fasta file')
    parser.add_argument('-l', '--readlength', dest='readlength', required=True, help='read length')
    parser.add_argument('--picard_path', dest='picard_path', required=True,
                        help='picard path directory')
    parser.add_argument('-o', '--outdir', dest='outdir', required=True, help='out directory(default current path)')

    parser.add_argument('--aligner', dest='aligner', required=True, default="bwa",
                        help='choose a aligner(default bwa)')
    parser.add_argument('--alignerIndex', dest='alignerIndex', required=True, default=0.9,
                        help='aligner genome index')
    parser.add_argument('--seqer', dest='seqer', required=True, help='seqer (illumina, life, BGI)')

    parser.add_argument('--floworder', dest='floworder', required=False, help='flower order of life sequence')
    parser.add_argument('--libkey', dest='libkey', required=False, help='libkey of life sequence')
    parser.add_argument('--barcode', dest='barcode', required=False, help='barcode of life sequence')
    parser.add_argument('--mindepth', dest='mindepth', required=False, default=30,
                        help='minimum depth(default = 1000)')
    parser.add_argument('-g', '--single', action='store_true', default=False,
                        help="input BAM is simgle-ended (default is paired-end)")
    parser.add_argument('-p', '--process', dest='process', required=False, default=1,
                        help='process num(default = 1)')
    parser.add_argument('--minmutreads', dest='minmutreads', required=False, default=5,
                        help='minimum depth(default = 5)')
    parser.add_argument('--tag', action='store_true', default=False, help="tag mutated reads (default False)")
    parser.add_argument('--minmapq', dest='minmapq', required=False, default=20,
                        help='read mapping quality less than minmapq will not be considered to edit')
    parser.add_argument('--multmapfilter', action='store_true', default=False,
                        help="multiple mapped reads will not be considered to edit. (default is True)")
    run_args = parser.parse_args()
    assert os.path.exists(run_args.svfile), "svfile is not existed."
    assert os.path.exists(run_args.bamfile), "bamfile is not existed."
    assert os.path.exists(run_args.bamfile + ".bai"), "index of bam file is not existed."
    assert os.path.exists(run_args.reffasta), "reference fasta file is not existed."
    if run_args.debug:
        print type(run_args)
    main(run_args)


if __name__ == "__main__":
    run()
