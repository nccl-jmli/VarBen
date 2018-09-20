from varben.common.bamobject import SV
import pysam
import traceback


def check_sv(sv_info, chrom_len):
    min_size = 100
    chr_id, start_id, end_id, type_id, freq_id = 0, 1, 2, 3, 4
    trans_chr_id, trans_start_id, trans_end_id = 5, 6, 7
    cnv_type_id, dup_num_id = 5, 5

    try:
        chrom, start, end = sv_info[chr_id], int(sv_info[start_id]), int(sv_info[end_id])
        sv_type = sv_info[type_id].lower()
        freq = float(sv_info[freq_id])
        if chrom not in chrom_len:
            raise Exception("chrom is not existed!")
        if not (1 <= start <= chrom_len[sv_info[chr_id]] and 1 <= end <= chrom_len[sv_info[chr_id]] and start <= end):
            raise Exception("start & end is not valid!")
        if not sv_type == "cnv" and not 0.0 < freq <= 1.0:
            raise Exception("freq is not valid!")
        if sv_type == "inv":
            if end - start < min_size:
                raise Exception("inversion is not in correct format")
        elif sv_type == "ins":
            pass
        elif sv_type == "del":
            if end - start < min_size:
                raise Exception("deletion is not in correct format")
        elif sv_type == "dup":
            print sv_info[dup_num_id]
            dup_num = int(sv_info[dup_num_id])
            if end - start < min_size or not dup_num >= 1:
                raise Exception("duplication is not in correct format")
        elif sv_type == "cnv":
            if (sv_info[cnv_type_id] == "loss" and float(sv_info[freq_id]) <= 1) or (
                    sv_info[cnv_type_id] == "gain" and float(sv_info[freq_id]) >= 1):
                pass
            else:
                raise Exception("cnv is not in correct format")
        elif sv_type == "trans_balance":
            trans_chr, trans_start, trans_end = sv_info[trans_chr_id], int(sv_info[trans_start_id]), int(
                sv_info[trans_end_id])
            if trans_chr not in chrom_len or end - start < 20 or trans_end - trans_start < 20:
                raise Exception("translocation_balance is not in correct format!")
        elif sv_type == "trans_unbalance":
            trans_chr, trans_start, trans_end = sv_info[trans_chr_id], int(sv_info[trans_start_id]), int(
                sv_info[trans_end_id])
            if trans_chr not in chrom_len or end - start < min_size or trans_end - trans_start != 1:
                raise Exception("translocation_unbalance is not in correct format!")
        elif sv_type == "trans_chrom":
            trans_chr, trans_start, trans_end = sv_info[trans_chr_id], int(sv_info[trans_start_id]), int(
                sv_info[trans_end_id])
            if trans_chr not in chrom_len or end != start or trans_end != trans_start:
                raise Exception("translocation_unbalance is not in correct format!")
        else:
            raise Exception("mutation type is not correct!")

        sv = SV(sv_info[chr_id], int(sv_info[start_id]) - 1, int(sv_info[end_id]) - 1, sv_type, float(sv_info[freq_id]))
        if sv_type in ("trans_balance", "trans_unbalance", "trans_chrom"):
            sv.trans_chrom = sv_info[trans_chr_id]
            sv.trans_start = int(sv_info[trans_start_id]) - 1
            sv.trans_end = int(sv_info[trans_end_id]) - 1
        elif sv_type == "cnv":
            sv.cnv_type = sv_info[cnv_type_id]
        elif sv_type == "dup":
            sv.dup_num = int(sv_info[dup_num_id])
        return True, sv
    except Exception as e:
        print traceback.format_exc()
        print e.args
        return False, e.args


def check_sv_file(sv_file, ref_file, invalid_log):
    ref = pysam.FastaFile(ref_file)
    chroms = ref.references
    chrom_len = {}
    for chrom in chroms:
        chrom_len[chrom] = ref.get_reference_length(chrom)
    sv_list = []
    with open(sv_file, 'r') as fin:
        for line in fin:
            line = line.rstrip()
            if line[0] == "#":
                continue
            if line == "":
                break
            sv_info = line.split()
            print sv_info
            flag, res = check_sv(sv_info, chrom_len)
            if not flag:
                invalid_log.info(line, res[0])
            sv_list.append(res)
    fin.close()
    return sv_list
