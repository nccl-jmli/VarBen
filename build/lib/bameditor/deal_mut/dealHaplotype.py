from bameditor.deal_mut.readEditor import editRead
from bameditor.common.methods import count_coverage, getReadStrand, find_mate
from bameditor.common.bamconvert import remap
from bameditor.common.bamobject import Read
from collections import OrderedDict
import random
import pysam
from multiprocessing import Pool
import os


def deal_haplotype(bam_file, haplotype, reffasta, haplotype_prefix, mindepth, minmutreads, minmapq, diffcover,
                   is_single, is_multmapfilter, aligner, aligner_index, invalid_log, **kwargs):
    reads_dict = OrderedDict()
    bam = pysam.AlignmentFile(bam_file, 'rb')
    reads = bam.fetch(reference=haplotype.chrom, start=haplotype.start, end=haplotype.end + 1)
    depth = 0
    for read in reads:
        depth += 1
        if read.reference_start is not None and not read.is_secondary and bin(read.flag & 2048) != bin(2048):
            if read.query_name not in reads_dict:
                reads_dict[read.query_name] = {}
            strand = getReadStrand(read)
            reads_dict[read.query_name][strand] = read

    # judge depth and mut reads whether qualified
    if depth < int(mindepth):
        print "depth less than min depth!"
        invalid_log.info("depth of haplosize in position %s:%s-%s les than min depth(%s)" % (
            haplotype.chrom, haplotype.start, haplotype.end, mindepth))
        return False
    else:
        mut_reads_num = int(depth * haplotype.freq)
        if mut_reads_num < int(minmutreads):
            print "mutation reads num less than minmutreads!"
            invalid_log.info("mut reads of haplosize in position %s:%s-%s less than min mut reads(%s)" % (
                haplotype.chrom, haplotype.start, haplotype.end, minmutreads))
            return False

    print "start pick reads"
    chosen_reads, mate_reads = pick_reads(bam, reads_dict, mut_reads_num, is_single, minmapq, is_multmapfilter)
    print "end pick reads"
    # edit
    my_chosen_reads = {}
    my_mate_reads = {}
    tmp_bam_file = haplotype_prefix + ".chosen.edited.bam"
    tmp_bam = pysam.AlignmentFile(tmp_bam_file, 'wb', template=bam)
    chosen_reads_num = 0

    real_mut_reads_num = 0
    for readName, readInfo in chosen_reads.items():
        my_chosen_reads[readName] = {}
        tmp_dict = {}
        tmp_dict2 = {}
        for strand, read in readInfo.items():
            my_read = Read(read)
            res = editRead(my_read, reffasta, haplotype.mutList)
            if res is False:
                continue
            real_mut_reads_num += 1
            sequence, quality, shift = res
            read.query_sequence = sequence
            read.query_qualities = quality
            tmp_dict[strand] = my_read
            tmp_dict2[strand] = read
        if is_single:
            for strand in tmp_dict:
                my_chosen_reads[readName][strand] = tmp_dict[strand]
                tmp_bam.write(tmp_dict2[strand])
                chosen_reads_num += 1
        else:
            if len(tmp_dict) == 0:
                continue
            elif len(tmp_dict) == 1 and readName in mate_reads:
                for strand in tmp_dict:
                    my_chosen_reads[readName][strand] = tmp_dict[strand]
                    tmp_bam.write(tmp_dict2[strand])
                    chosen_reads_num += 1
                mate_read = mate_reads[readName]
                my_mate_reads[readName] = Read(mate_read)
                tmp_bam.write(mate_read)
            elif len(tmp_dict) == 2:
                for strand in tmp_dict:
                    my_chosen_reads[readName][strand] = tmp_dict[strand]
                    tmp_bam.write(tmp_dict2[strand])
                    chosen_reads_num += 1
    tmp_bam.close()

    # alignment and judge coverdiff whether qualified
    chosen_bam_file = haplotype_prefix + ".chosen.remap.bam"
    genome_index = aligner_index
    remap(genome_index, tmp_bam_file, chosen_bam_file, aligner, is_single)
    chosen_bam = pysam.AlignmentFile(chosen_bam_file)
    if judge_coverdiff(bam, depth, chosen_bam, chosen_reads_num, haplotype, float(diffcover)):
        return my_chosen_reads, my_mate_reads, real_mut_reads_num, depth
    else:
        print "coverdiff is less than minDiffCover"
        return False


def pick_reads(bam, reads_dict, choose_num, is_single, minmapq, is_multmapfilter):
    total_reads_num = len(reads_dict)
    print total_reads_num, choose_num
    chosen_reads = {}
    keys = reads_dict.keys()
    chosen_reads_id = []
    mate_reads = {}
    if is_single:
        i = 0
        chosen_reads_id = []
        while i < choose_num:
            choose_id = random.randint(0, total_reads_num - 1)
            if choose_id in chosen_reads_id:
                continue
            read_name = keys[choose_id]
            strand_keys = reads_dict[read_name].keys()
            strand = strand_keys[0]
            read = reads_dict[read_name][strand]
            if int(read.mapping_quality) < minmapq or len(read.query_sequence) < 250:
                continue
            chosen_reads[read_name] = {}
            chosen_reads[read_name][strand] = read
            chosen_reads_id.append(choose_id)
            i += 1

    else:
        num = 0
        while num < choose_num:
            choose_id = random.randint(0, total_reads_num - 1)
            if choose_id in chosen_reads_id:
                continue
            read_name = keys[choose_id]

            strand_keys = reads_dict[read_name].keys()
            if len(reads_dict[read_name]) == 1:
                read = reads_dict[read_name][strand_keys[0]]
                if int(read.mapping_quality) < minmapq:
                    continue
                else:
                    if is_multmapfilter and read.has_tag("XA"):
                        continue
                    chosen_reads[read_name] = {}
                    chosen_reads[read_name][strand_keys[0]] = read
                    if not read.mate_is_unmapped:
                        mate_read = find_mate(read, bam)
                        if mate_read is not None:
                            mate_reads[read_name] = mate_read
                        else:
                            print "no mate read found!"
            if len(reads_dict[read_name]) == 2:
                read1 = reads_dict[read_name][strand_keys[0]]
                read2 = reads_dict[read_name][strand_keys[1]]
                if int(read1.mapping_quality) < minmapq or int(read2.mapping_quality) < minmapq:
                    continue
                else:
                    if is_multmapfilter and (read1.has_tag("XA") or read2.has_tag("XA")):
                        continue
                    chosen_reads[read_name] = {}
                    chosen_reads[read_name][strand_keys[0]] = read1
                    chosen_reads[read_name][strand_keys[1]] = read2

            num += len(reads_dict[read_name])
            chosen_reads_id.append(choose_id)
    return chosen_reads, mate_reads


def judge_coverdiff(bam, depth, chosen_bam, chosen_reads_num, haplotype, min_diff_cover):
    old_coverage = count_coverage(bam, haplotype.chrom, haplotype.start, haplotype.end)
    new_coverage = count_coverage(chosen_bam, haplotype.chrom, haplotype.start, haplotype.end) + (
        depth - chosen_reads_num)
    if new_coverage * 1.0 / old_coverage > min_diff_cover:
        return True
    else:
        return False


def deal_haplotype_multi(bam_file, haplotype_list, out_dir, reffasta, process, mindepth,
                         minmutreads, minmapq, diffcover, is_single, is_multmapfilter,
                         aligner, aligner_index, invalid_log):
    haplotype_pool = Pool(processes=int(process))
    haplotype_res = []
    haplotype_temp_out_dir = os.path.join(out_dir, "haplotype_out")

    for haplotype in haplotype_list:
        haplotype_prefix = os.path.join(haplotype_temp_out_dir,
                                        "%s_%s_%s" % (haplotype.chrom, haplotype.start, haplotype.end))
        haplotype_res.append(
            haplotype_pool.apply_async(deal_haplotype,
                                       args=(bam_file, haplotype, reffasta, haplotype_prefix, mindepth,
                                             minmutreads, minmapq, diffcover, is_single, is_multmapfilter,
                                             aligner, aligner_index, invalid_log)))
        haplotype_pool.close()
        haplotype_pool.join()

    # step3: merge mut_list of each read
    total_chosen_reads = {}
    total_chosen_reads_muts = {}
    total_chosen_mate_reads = {}
    haplotype_out_file = os.path.join(out_dir, "haplotype.txt")
    hap_out = open(haplotype_out_file, 'w')

    for each_res, haplotype in zip(haplotype_res, haplotype_list):
        res = each_res.get()
        if not res:
            continue
        chosen_reads, mate_reads, real_mut_reads_num, depth = res
        hap_out.write(
            "\t".join(
                [str(haplotype), str(depth), str(real_mut_reads_num),
                 str(real_mut_reads_num * 1.0 / depth)]) + "\n")
        # invalid_log.info(invalid)
        for read_name in chosen_reads:
            if read_name not in total_chosen_reads:
                total_chosen_reads[read_name] = {}
                total_chosen_reads_muts[read_name] = {}
            for strand in chosen_reads[read_name]:
                if strand not in total_chosen_reads[read_name]:
                    total_chosen_reads[read_name][strand] = chosen_reads[read_name][strand]
                    total_chosen_reads_muts[read_name][strand] = []
                total_chosen_reads_muts[read_name][strand].extend(haplotype.mutList)

        for read_name in mate_reads:
            if read_name not in total_chosen_mate_reads:
                total_chosen_mate_reads[read_name] = mate_reads[read_name]
        hap_out.close()

        return total_chosen_reads, total_chosen_reads_muts
