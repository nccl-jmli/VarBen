from bameditor.common.methods import count_coverage
from bameditor.deal_sv.dealReadsType import deal_type1, deal_type2, deal_type3, deal_type4, deal_type5, deal_type6, \
    deal_type7
from bameditor.deal_sv.getReadsType import pos_type_classify
from bameditor.common.methods import getKeyName
import pysam
import os


def get_write_reads(total_modify_reads, total_delete_reads, total_add_reads, total_reads_file_dict,
                    total_reads_list_dict):
    for typ, reads_dict in zip(('modify', 'delete', 'add'), (total_modify_reads, total_delete_reads, total_add_reads)):
        reads_file = total_reads_file_dict[typ]
        reads_pair = total_reads_list_dict[typ]
        for read_pair in reads_dict:
            tmp = []
            for read in read_pair:
                read_name = getKeyName(read)
                tmp.append(read_name)
                reads_file.write(str(read))
            reads_pair.write("%s\n" % ",".join(tmp))


def deal_sv(bam_file, ref_file, sv_list, is_single, minmapq, is_multmapfilter, mindepth, minmutreads, read_length,
            out_dir, insert_size,
            invalid_log, run_log, success_list):
    ref = pysam.FastaFile(ref_file)
    # total_modify_reads, total_delete_reads, total_add_reads = [], [], []
    bam = pysam.AlignmentFile(bam_file, 'rb')
    total_reads_file_dict = {}
    total_reads_list_dict = {}
    reads_pair_list, reads_file_list = [], []
    fout_success = open(success_list, 'w')
    for typ in ('modify', 'delete', 'add'):
        reads_file = os.path.join(out_dir, "total_%s_reads.txt" % typ)
        total_reads_file_dict[typ] = open(reads_file, 'w')
        pair_list = os.path.join(out_dir, "total_%s_pair.list" % typ)
        total_reads_list_dict[typ] = open(pair_list, 'w')
        reads_file_list.append(reads_file)
        reads_pair_list.append(pair_list)

    for sv_info in sv_list:
        freq = sv_info.freq
        chrom = sv_info.chrom
        start = sv_info.start
        end = sv_info.end
        mindepth = int(mindepth)
        minmutreads = int(minmutreads)
        if not sv_info.sv_type == "cnv":
            start_coverage = count_coverage(bam, chrom, start, start)
            start_coverage_deal = start_coverage * freq
            end_coverage = count_coverage(bam, chrom, end, end)
            end_coverage_deal = end_coverage
            if start_coverage < mindepth or end_coverage < mindepth:
                invalid_log.info(str(sv_info), "coverage too small")
                continue
            if start_coverage_deal < minmutreads or end_coverage_deal < minmutreads:
                invalid_log.info(str(sv_info), "coverage deal too small")
                continue

        func_dict = {"del": deal_deletion, "inv": deal_inversion, "ins": deal_insertion,
                     "trans_balance": deal_translocation_balance,
                     "trans_unbalance": deal_translocation_unbalance, "trans_chrom": deal_translocation_chrom,
                     "dup": deal_duplication, "cnv": deal_cnv}
        func = func_dict[sv_info.sv_type]
        run_log.info("=====> " + str(sv_info))
        res = func(bam_file, ref, sv_info, is_single, read_length, out_dir, insert_size, minmapq, is_multmapfilter)
        if res[0] is False:
            invalid_log.info(str(sv_info), res[1])
            continue
        fout_success.write(sv_info.svinfo()+"\n")
        modify_reads, delete_reads, add_reads = res
        # total_modify_reads.extend(modify_reads)
        # total_delete_reads.extend(delete_reads)
        # total_add_reads.extend(add_reads)
        get_write_reads(modify_reads, delete_reads, add_reads, total_reads_file_dict,
                        total_reads_list_dict)
    fout_success.close()
    bam.close()

    for typ in ('modify', 'delete', 'add'):
        total_reads_file_dict[typ].close()
        total_reads_list_dict[typ].close()

    return reads_file_list + reads_pair_list
    # return total_modify_reads, total_delete_reads, total_add_reads


def deal_insertion(bam, ref, svinfo, is_single, readLength, tempDir, insertSize, minmapq, is_multmapfilter):
    pass


def deal_translocation_unbalance(bam, ref, svinfo, is_single, readLength, tempDir, insertSize, minmapq,
                                 is_multmapfilter):
    freq = svinfo.freq
    chrom = svinfo.chrom
    start = svinfo.start
    end = svinfo.end
    trans_chr, trans_start, trans_end = svinfo.trans_chrom, svinfo.trans_start, svinfo.trans_end
    info = "%s:%s-%s" % (chrom, start, end)
    svtype = "trans_unbalance"
    if is_single:
        reads_type6_left, reads_type6_right, reads_type7, filtered_read_num = pos_type_classify(bam, chrom, start, end,
                                                                                                is_single, readLength,
                                                                                                tempDir)
        trans_reads_type6_left, trans_reads_type6_right, trans_reads_type7, filtered_read_num = pos_type_classify(bam,
                                                                                                                  trans_chr,
                                                                                                                  trans_start,
                                                                                                                  trans_end,
                                                                                                                  is_single,
                                                                                                                  readLength,
                                                                                                                  tempDir,
                                                                                                                  minmapq=minmapq,
                                                                                                                  is_multmapfilter=is_multmapfilter)
        type_reads_num = len(reads_type6_left) + len(reads_type6_right) + len(reads_type7)
        freq = freq * (type_reads_num + filtered_read_num) / type_reads_num
        if freq > 1:
            freq = 1
        modify_reads_type6 = deal_type6(ref, trans_reads_type6_left, trans_reads_type6_left, freq, trans_start,
                                        trans_end, svtype, subPos=info)
        add_reads_type7 = deal_type7(reads_type7, freq, svtype)

        return modify_reads_type6, [], add_reads_type7

    else:
        reads_type1_left, reads_type1_right, reads_type2_left, reads_type2_right, reads_type3_left, reads_type3_right, reads_type4, reads_type5_left, reads_type5_right, filtered_read_num = pos_type_classify(
            bam, chrom, start, end, is_single, readLength, tempDir)
        trans_reads_type1_left, trans_reads_type1_right, trans_reads_type2_left, trans_reads_type2_right, trans_reads_type3_left, trans_reads_type3_right, trans_reads_type4, trans_reads_type5_left, trans_reads_type5_right, filtered_read_num = pos_type_classify(
            bam, trans_chr, trans_start, trans_end, is_single, readLength, tempDir, extension=insertSize[1],
            minmapq=minmapq, is_multmapfilter=is_multmapfilter)
        # exit()
        type_reads_num = len(reads_type1_left) + len(reads_type1_right) + len(reads_type2_left) + len(
            reads_type2_right) + len(reads_type3_left) + len(reads_type3_right) + len(reads_type4) + len(
            reads_type5_left) + len(reads_type5_right)
        if type_reads_num == 0:
            return False, "All reads is filtered in this scope, total & filtered reads num %s" % filtered_read_num
        freq = freq * (type_reads_num * 2 + filtered_read_num) / (type_reads_num * 2)

        modify_reads_type1 = deal_type1(ref, trans_reads_type1_left, [], freq, trans_start, trans_end, svtype,
                                        subPos=info)
        modify_reads_type2 = deal_type2(ref, [], trans_reads_type2_left, freq, trans_start, trans_end, svtype,
                                        subPos=info)
        res = deal_type3(trans_reads_type3_left,
                         trans_reads_type3_left, freq, insertSize,
                         trans_start,
                         trans_end, svtype,
                         supple1=reads_type3_left,
                         supple2=reads_type3_right,
                         subPos=info)
        if res is False:
            return False, "The reads are not satisfied to modify, see more details in log file"
        modify_reads_type3, add_reads_type3, delete_reads_type3 = res
        add_reads_type4 = deal_type4(reads_type4, freq, svtype, insertSize=insertSize)
        modify_reads_type5 = deal_type5(ref, trans_reads_type5_left, trans_reads_type5_left, freq, trans_start,
                                        trans_end, svtype, subPos=info)

        return modify_reads_type1 + modify_reads_type2 + modify_reads_type3 + modify_reads_type5, [], add_reads_type4


def deal_cnv(bam, ref, svinfo, is_single, readLength, tempDir, insertSize, minmapq, is_multmapfilter):
    freq = svinfo.freq
    chrom = svinfo.chrom
    start = svinfo.start
    end = svinfo.end
    svtype = "cnv"
    if is_single:
        reads_type6_left, reads_type6_right, reads_type7, filtered_read_num = pos_type_classify(bam, chrom, start, end,
                                                                                                is_single, readLength,
                                                                                                tempDir,
                                                                                                minmapq=minmapq,
                                                                                                is_multmapfilter=is_multmapfilter)
        type_reads_num = len(reads_type6_left) + len(reads_type6_right) + len(reads_type7)
        freq = freq * (type_reads_num + filtered_read_num) / type_reads_num
        if svinfo.cnv_type == "loss":
            delete_reads_type7 = deal_type7(reads_type7, freq, svtype)
            return [], delete_reads_type7, []
        elif svinfo.cnv_type == "gain":
            add_reads_type7 = deal_type7(reads_type7, freq - 1, svtype)
            return [], [], add_reads_type7

    else:
        reads_type1_left, reads_type1_right, reads_type2_left, reads_type2_right, reads_type3_left, reads_type3_right, reads_type4, reads_type5_left, reads_type5_right, filtered_read_num = pos_type_classify(
            bam, chrom, start, end, is_single, readLength, tempDir, minmapq=minmapq, is_multmapfilter=is_multmapfilter)
        type_reads_num = len(reads_type1_left) + len(reads_type1_right) + len(reads_type2_left) + len(
            reads_type2_right) + len(reads_type3_left) + len(reads_type3_right) + len(reads_type4) + len(
            reads_type5_left) + len(reads_type5_right)
        if type_reads_num == 0:
            return False, "All reads is filtered in this scope, total & filtered reads num %s" % filtered_read_num
        if svinfo.cnv_type == "loss":
            freq = freq * (type_reads_num * 2 + filtered_read_num) / (type_reads_num * 2)
            delete_reads_type4 = deal_type4(reads_type4, freq, svtype, insertSize=insertSize, cnvType="loss")
            return [], delete_reads_type4, []
        elif svinfo.cnv_type == "gain":
            freq = (freq - 1) * (type_reads_num * 2 + filtered_read_num) / (type_reads_num * 2)
            add_reads_type4 = deal_type4(reads_type4, freq, svtype, insertSize=insertSize, cnvType="gain")
            return [], [], add_reads_type4


def deal_duplication(bam, ref, svinfo, is_single, readLength, tempDir, insertSize, minmapq, is_multmapfilter):
    freq = svinfo.freq
    chrom = svinfo.chrom
    start = svinfo.start
    end = svinfo.end
    dup_num = svinfo.dup_num
    svtype = "dup"
    if is_single:
        reads_type6_left, reads_type6_right, reads_type7, filtered_read_num = pos_type_classify(bam, chrom, start, end,
                                                                                                is_single, readLength,
                                                                                                tempDir,
                                                                                                minmapq=minmapq,
                                                                                                is_multmapfilter=is_multmapfilter)
        type_reads_num = len(reads_type6_left) + len(reads_type6_right) + len(reads_type7)
        if type_reads_num == 0:
            return False, "All reads is filtered in this scope, total & filtered reads num %s" % filtered_read_num
        freq = freq * (type_reads_num + filtered_read_num) / type_reads_num
        if freq > 1:
            used_freq = 1
        else:
            used_freq = freq
        add_reads_type6 = deal_type6(ref, reads_type6_left, reads_type6_left, used_freq, start, end, svtype)
        add_reads_type7 = deal_type7(reads_type7, freq, svtype)
        return [], [], add_reads_type6 + add_reads_type7

    else:

        reads_type1_left, reads_type1_right, reads_type2_left, reads_type2_right, reads_type3_left, reads_type3_right, reads_type4, reads_type5_left, reads_type5_right, filtered_read_num = pos_type_classify(
            bam, chrom, start, end, is_single, readLength, tempDir, minmapq=minmapq, is_multmapfilter=is_multmapfilter)

        type_reads_num = len(reads_type1_left) + len(reads_type1_right) + len(reads_type2_left) + len(
            reads_type2_right) + len(reads_type3_left) + len(reads_type3_right) + len(reads_type4) + len(
            reads_type5_left) + len(reads_type5_right)
        if type_reads_num == 0:
            return False, "All reads is filtered in this scope, total & filtered reads num %s" % filtered_read_num
        freq = freq * (type_reads_num * 2 + filtered_read_num) / (type_reads_num * 2)
        print freq
        if freq > 1:
            freq = 1
        freq_dup = (dup_num - 1) * freq

        freq_dup_part = freq_dup / 2.0
        if freq_dup_part > 1:
            freq_dup_part = 1

        if freq_dup > 1:
            used_freq = 1
        else:
            used_freq = freq_dup
        # start_coverage = count_coverage(bam, chrom, start, start)
        # end_coverage = count_coverage(bam, chrom, end, end)
        add_reads = []

        add_reads_type2 = deal_type2(ref, reads_type2_left, reads_type2_right, used_freq, start, end, svtype)
        add_reads.extend(add_reads_type2)
        # modify by fangshs 20180606
        supple1 = reads_type3_right + reads_type1_right
        supple2 = reads_type3_left + reads_type1_left
        res = deal_type3(reads_type3_left, reads_type3_right,
                         freq_dup, insertSize, start, end, svtype, supple1=supple1, supple2=supple2)

        if res is False:
            return False, "The reads are not satisfied to modify, see more details in log file"
        modify_reads_type3, add_reads_type3, delete_reads_type3 = res
        add_reads.extend(add_reads_type3)
        add_reads_type4 = deal_type4(reads_type4, freq_dup, svtype, insertSize=insertSize)
        add_reads.extend(add_reads_type4)
        add_reads_type5 = deal_type5(ref, reads_type5_left, reads_type5_right, freq_dup_part, start, end, svtype)
        add_reads.extend(add_reads_type5)

        return [], [], add_reads


def deal_deletion(bam, ref, svinfo, is_single, readLength, tempDir, insertSize, minmapq, is_multmapfilter):
    freq = svinfo.freq
    chrom = svinfo.chrom
    start = svinfo.start
    end = svinfo.end
    svtype = "del"
    if is_single:
        reads_type6_left, reads_type6_right, reads_type7, filtered_read_num = pos_type_classify(bam, chrom, start, end,
                                                                                                is_single, readLength,
                                                                                                tempDir)
        type_reads_num = len(reads_type6_left) + len(reads_type6_right) + len(reads_type7)
        freq = freq * (type_reads_num + filtered_read_num) / type_reads_num
        if freq > 1:
            freq = 1
        modify_reads_type6 = deal_type6(ref, reads_type6_left, reads_type6_left, freq, start, end, svtype)
        delete_reads_type7 = deal_type7(reads_type7, freq, svtype)
        return modify_reads_type6, delete_reads_type7, []
    else:
        reads_type1_left, reads_type1_right, reads_type2_left, reads_type2_right, reads_type3_left, reads_type3_right, reads_type4, reads_type5_left, reads_type5_right, filtered_read_num = pos_type_classify(
            bam, chrom, start, end, is_single, readLength, tempDir, minmapq=minmapq, is_multmapfilter=is_multmapfilter)
        type_reads_num = len(reads_type1_left) + len(reads_type1_right) + len(reads_type2_left) + len(
            reads_type2_right) + len(reads_type3_left) + len(reads_type3_right) + len(reads_type4) + len(
            reads_type5_left) + len(reads_type5_right)
        if type_reads_num == 0:
            return False, "All reads is filtered in this scope, total & filtered reads num %s" % filtered_read_num
        if type_reads_num == 0:
            return False, "All reads is filtered in this scope, total & filtered reads num %s" % filtered_read_num
        freq = freq * (type_reads_num * 2 + filtered_read_num) / (type_reads_num * 2)
        if freq > 1:
            freq = 1

        modify_reads_type1 = deal_type1(ref, reads_type1_left, reads_type1_right, freq, start, end, svtype)
        delete_reads_type2 = deal_type2(ref, reads_type2_left, reads_type2_right, freq, start, end, svtype)
        supple1 = reads_type3_right + reads_type1_right
        res = deal_type3(reads_type3_left, reads_type3_right, freq, insertSize, start, end, svtype, supple1=supple1)
        if res is False:
            return False, "The reads are not satisfied to modify, see more details in log file"

        modify_reads_type3, add_reads_type3, delete_reads_type3 = res
        delete_reads_type4 = deal_type4(reads_type4, freq, svtype)
        delete_reads_type5 = deal_type5(ref, reads_type5_left, reads_type5_right, freq, start, end, svtype)

        # modify, delete, add
        return modify_reads_type1 + modify_reads_type3, delete_reads_type2 + delete_reads_type4 + delete_reads_type5 + delete_reads_type3, []


def deal_inversion(bam, ref, svinfo, is_single, readLength, tempDir, insertSize, minmapq, is_multmapfilter):
    freq = svinfo.freq
    chrom = svinfo.chrom
    start = svinfo.start
    end = svinfo.end
    svtype = "inv"
    if is_single:
        reads_type6_left, reads_type6_right, reads_type7, filtered_read_num = pos_type_classify(bam, chrom, start, end,
                                                                                                is_single, readLength,
                                                                                                tempDir, center=False,
                                                                                                maxsize=1,
                                                                                                minmapq=minmapq,
                                                                                                is_multmapfilter=is_multmapfilter)
        type_reads_num = len(reads_type6_left) + len(reads_type6_right) + len(reads_type7)
        if type_reads_num == 0:
            return False, "All reads is filtered in this scope, total & filtered reads num %s" % filtered_read_num
        freq = freq * (type_reads_num + filtered_read_num) / type_reads_num
        if freq > 1:
            freq = 1
        modify_reads_type6 = deal_type6(ref, reads_type6_left, reads_type6_left, freq, start, end, svtype)
        return modify_reads_type6, [], []
    else:
        reads_type1_left, reads_type1_right, reads_type2_left, reads_type2_right, reads_type3_left, reads_type3_right, reads_type4, reads_type5_left, reads_type5_right, filtered_read_num = pos_type_classify(
            bam, chrom, start, end, is_single, readLength, tempDir, center=False, maxsize=insertSize[1],
            minmapq=minmapq, is_multmapfilter=is_multmapfilter)
        type_reads_num = len(reads_type1_left) + len(reads_type1_right) + len(reads_type2_left) + len(
            reads_type2_right) + len(reads_type3_left) + len(reads_type3_right) + len(reads_type4) + len(
            reads_type5_left) + len(reads_type5_right)
        if type_reads_num == 0:
            return False, "All reads is filtered in this scope, total & filtered reads num %s" % filtered_read_num
        freq = freq * (type_reads_num * 2 + filtered_read_num) / (type_reads_num * 2)
        if freq > 1:
            freq = 1
        modified_reads = []
        modify_reads_type1 = deal_type1(ref, reads_type1_left, reads_type1_right, freq, start, end, svtype)
        modified_reads.extend(modify_reads_type1)
        modify_reads_type2 = deal_type2(ref, reads_type2_left, reads_type2_right, freq, start, end, svtype)
        modified_reads.extend(modify_reads_type2)
        res = deal_type3(reads_type3_left, reads_type3_right, freq, insertSize, start, end, svtype, supple1=reads_type4)
        if res is False:
            return False, "The reads are not satisfied to modify, see more details in log file."

        modify_reads_type3, add_reads_type3, delete_reads_type3 = res
        modified_reads.extend(modify_reads_type3)
        modify_reads_type5 = deal_type5(ref, reads_type5_left, reads_type5_right, freq, start, end, svtype)
        modified_reads.extend(modify_reads_type5)

        return modified_reads, [], []


def deal_translocation_balance(bam, ref, svinfo, is_single, readLength, tempDir, insertSize, minmapq, is_multmapfilter):
    freq = svinfo.freq
    chrom = svinfo.chrom
    start = svinfo.start
    end = svinfo.end
    trans_chr, trans_start, trans_end = svinfo.trans_chrom, svinfo.trans_start, svinfo.trans_end
    trans_info = "%s:%s-%s" % (trans_chr, trans_start, trans_end)
    info = "%s:%s-%s" % (chrom, start, end)
    svtype = "trans_balance"
    if is_single:
        reads_type6_left, reads_type6_right, reads_type7, filtered_read_num = pos_type_classify(bam, chrom, start, end,
                                                                                                is_single, readLength,
                                                                                                tempDir, center=False,
                                                                                                maxsize=1,
                                                                                                minmapq=minmapq,
                                                                                                is_multmapfilter=is_multmapfilter)
        trans_reads_type6_left, trans_reads_type6_right, trans_reads_type7, filtered_read_num = pos_type_classify(bam,
                                                                                                                  trans_chr,
                                                                                                                  trans_start,
                                                                                                                  trans_end,
                                                                                                                  is_single,
                                                                                                                  readLength,
                                                                                                                  tempDir,
                                                                                                                  minmapq=minmapq,
                                                                                                                  is_multmapfilter=is_multmapfilter)
        type_reads_num = len(reads_type6_left) + len(reads_type6_right) + len(reads_type7)
        if type_reads_num == 0:
            return False, "All reads is filtered in this scope, total & filtered reads num %s" % filtered_read_num
        freq = freq * (type_reads_num + filtered_read_num) / type_reads_num
        if freq > 1:
            freq = 1
        modify_reads_type6 = deal_type6(ref, reads_type6_left, reads_type6_right, freq, start, end, svtype,
                                        subPos=trans_info)
        trans_modify_reads_type6 = deal_type6(ref, trans_reads_type6_left, trans_reads_type6_right, freq, trans_start,
                                              trans_end, svtype, subPos=info)
        return modify_reads_type6 + trans_modify_reads_type6, [], []

    else:
        reads_type1_left, reads_type1_right, reads_type2_left, reads_type2_right, reads_type3_left, reads_type3_right, reads_type4, reads_type5_left, reads_type5_right, filtered_read_num = pos_type_classify(
            bam, chrom, start, end, is_single, readLength, tempDir + ".info", center=False, maxsize=insertSize[1],
            minmapq=minmapq, is_multmapfilter=is_multmapfilter)
        # exit()
        trans_reads_type1_left, trans_reads_type1_right, trans_reads_type2_left, trans_reads_type2_right, trans_reads_type3_left, trans_reads_type3_right, trans_reads_type4, trans_reads_type5_left, trans_reads_type5_right, filtered_read_num = pos_type_classify(
            bam, trans_chr, trans_start, trans_end, is_single, readLength, tempDir + ".trans",
            center=False, maxsize=insertSize[1], minmapq=minmapq, is_multmapfilter=is_multmapfilter)
        type_reads_num = len(reads_type1_left) + len(reads_type1_right) + len(reads_type2_left) + len(
            reads_type2_right) + len(reads_type3_left) + len(reads_type3_right) + len(reads_type4) + len(
            reads_type5_left) + len(reads_type5_right)
        if type_reads_num == 0:
            return False, "All reads is filtered in this scope, total & filtered reads num %s" % filtered_read_num
        freq = freq * (type_reads_num * 2 + filtered_read_num) / (type_reads_num * 2)
        if freq > 1:
            freq = 1
        modified_reads = []
        modify_reads_type1 = deal_type1(ref, reads_type1_left, reads_type1_right, freq, start, end, svtype,
                                        subPos=trans_info)
        trans_modify_reads_type1 = deal_type1(ref, trans_reads_type1_left, trans_reads_type1_right, freq, trans_start,
                                              trans_end, svtype, subPos=info)
        modified_reads.extend(modify_reads_type1 + trans_modify_reads_type1)
        modify_reads_type2 = deal_type2(ref, reads_type2_left, reads_type2_right, freq, start, end, svtype,
                                        subPos=trans_info)
        trans_modify_reads_type2 = deal_type2(ref, trans_reads_type2_left, trans_reads_type2_right, freq, trans_start,
                                              trans_end, svtype, subPos=info)
        modified_reads.extend(modify_reads_type2 + trans_modify_reads_type2)
        res = deal_type3(reads_type3_left, reads_type3_right, freq,
                         insertSize, start, end, svtype,
                         subPos=trans_info,
                         supple1=trans_reads_type3_left,
                         supple2=trans_reads_type3_right)
        res_trans = deal_type3(trans_reads_type3_left, trans_reads_type3_right,
                               freq, insertSize,
                               trans_start, trans_end,
                               svtype, subPos=info,
                               supple1=reads_type3_left,
                               supple2=reads_type3_right)
        if res is False or res_trans is False:
            return False, "The reads are not satisfied to modify, see more details in log file"

        modify_reads_type3, add_reads_type3, delete_reads_type3 = res
        trans_modify_reads_type3, trans_add_reads_type3, trans_delete_reads_type3 = res_trans
        modified_reads.extend(modify_reads_type3 + trans_modify_reads_type3)
        modify_reads_type5 = deal_type5(ref, reads_type5_left, reads_type5_right, freq, start, end, svtype,
                                        subPos=trans_info)
        trans_modify_reads_type5 = deal_type5(ref, trans_reads_type5_left, trans_reads_type5_right, freq, trans_start,
                                              trans_end, svtype, subPos=info)
        modified_reads.extend(modify_reads_type5 + trans_modify_reads_type5)

        return modified_reads, [], []


def deal_translocation_chrom(bam, ref, svinfo, is_single, read_length, tempDir, insertSize, minmapq, is_multmapfilter):
    freq = svinfo.freq
    chrom = svinfo.chrom
    start = svinfo.start
    end = svinfo.end
    trans_chr, trans_start, trans_end = svinfo.trans_chrom, svinfo.trans_start, svinfo.trans_end
    trans_info = "%s:%s-%s" % (trans_chr, trans_start, trans_end)
    info = "%s:%s-%s" % (chrom, start, end)
    svtype = "trans_chrom"
    if is_single:
        pass

    else:
        reads_type1_left, reads_type1_right, reads_type2_left, reads_type2_right, reads_type3_left, reads_type3_right, reads_type4, reads_type5_left, reads_type5_right, filtered_read_num = pos_type_classify(
            bam, chrom, start, end, is_single, read_length, tempDir + ".info", center=False, maxsize=insertSize[1],
            minmapq=minmapq, is_multmapfilter=is_multmapfilter)
        # exit()
        trans_reads_type1_left, trans_reads_type1_right, trans_reads_type2_left, trans_reads_type2_right, trans_reads_type3_left, trans_reads_type3_right, trans_reads_type4, trans_reads_type5_left, trans_reads_type5_right, filtered_read_num = pos_type_classify(
            bam, trans_chr, trans_start, trans_end, is_single, read_length, tempDir + ".trans",
            center=False, maxsize=insertSize[1], minmapq=minmapq, is_multmapfilter=is_multmapfilter)
        type_reads_num = len(reads_type1_left) + len(reads_type1_right) + len(reads_type2_left) + len(
            reads_type2_right) + len(reads_type3_left) + len(reads_type3_right) + len(reads_type4) + len(
            reads_type5_left) + len(reads_type5_right)
        if type_reads_num == 0:
            return False, "All reads is filtered in this scope, total & filtered reads num %s" % filtered_read_num
        freq = freq * (type_reads_num * 2 + filtered_read_num) / (type_reads_num * 2)
        if freq > 1:
            freq = 1
        modified_reads = []
        modify_reads_type1 = deal_type1(ref, reads_type1_left, reads_type1_right, freq, start, end, svtype,
                                        subPos=trans_info)
        trans_modify_reads_type1 = deal_type1(ref, trans_reads_type1_left, trans_reads_type1_right, freq, trans_start,
                                              trans_end, svtype, subPos=info)
        modified_reads.extend(modify_reads_type1 + trans_modify_reads_type1)
        modify_reads_type2 = deal_type2(ref, reads_type2_left, reads_type2_right, freq, start, end, svtype,
                                        subPos=trans_info)
        trans_modify_reads_type2 = deal_type2(ref, trans_reads_type2_left, trans_reads_type2_right, freq, trans_start,
                                              trans_end, svtype, subPos=info)
        modified_reads.extend(modify_reads_type2 + trans_modify_reads_type2)
        res = deal_type3(reads_type3_left, reads_type3_right, freq,
                         insertSize, start, end, svtype,
                         subPos=trans_info,
                         supple1=trans_reads_type3_left,
                         supple2=trans_reads_type3_right)
        res_trans = deal_type3(trans_reads_type3_left,
                               trans_reads_type3_right,
                               freq, insertSize,
                               trans_start, trans_end,
                               svtype, subPos=info,
                               supple1=reads_type3_left,
                               supple2=reads_type3_right)
        if res is False or res_trans is False:
            return False, "The reads are not satisfied to modify, see more details in log file"
        modify_reads_type3, add_reads_type3, delete_reads_type3 = res
        trans_modify_reads_type3, trans_add_reads_type3, trans_delete_reads_type3 = res_trans
        modified_reads.extend(modify_reads_type3 + trans_modify_reads_type3)
        modify_reads_type5 = deal_type5(ref, reads_type5_left, reads_type5_right, freq, start, end, svtype,
                                        subPos=trans_info)
        trans_modify_reads_type5 = deal_type5(ref, trans_reads_type5_left, trans_reads_type5_right, freq, trans_start,
                                              trans_end, svtype, subPos=info)
        modified_reads.extend(modify_reads_type5 + trans_modify_reads_type5)
        return modified_reads, [], []
