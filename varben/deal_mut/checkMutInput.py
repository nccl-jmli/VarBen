from varben.common.bamobject import Mutation, Haplotype
from varben.common.methods import calPolymorphism
import pysam


def check_alt(alt):
    for x in alt:
        if x not in 'ATCG':
            return False
    return True


def check_mut_info(mutinfo, ref):
    chroms = ref.references
    chrom_len = {}
    for chrom in chroms:
        chrom_len[chrom] = ref.get_reference_length(chrom)
    chr_id, start_id, end_id, vaf_id, type_id, alt_id = (ind for ind in range(6))
    try:
        if len(mutinfo) < 6:
            raise Exception("has less than 6 colums!")
        if mutinfo[chr_id] not in chrom_len:
            raise Exception("chrom is not existed!")
        mutinfo[start_id], mutinfo[end_id] = int(mutinfo[start_id]), int(mutinfo[end_id])
        if not (1 <= mutinfo[start_id] <= chrom_len[mutinfo[chr_id]] and 1 <= mutinfo[end_id] <= chrom_len[
            mutinfo[chr_id]] and mutinfo[start_id] <= chrom_len[mutinfo[chr_id]]):
            raise Exception("start & end is not valid!")
        mutinfo[vaf_id] = float(mutinfo[vaf_id])
        if not 0.001 <= mutinfo[vaf_id] <= 1:
            raise Exception("mutate frequency is not valid!")
        mut_type = mutinfo[type_id].lower()
        mutinfo[type_id] = mut_type
        if mut_type == "snv":
            if mutinfo[start_id] != mutinfo[end_id] or len(mutinfo[alt_id]) != 1 or not check_alt(mutinfo[alt_id]):
                raise Exception("not in %s correct format!" % mut_type)
            if ref.fetch(mutinfo[chr_id], mutinfo[start_id] - 1, mutinfo[end_id]) == mutinfo[alt_id]:
                raise Exception("alt base is equal to reference!")
        elif mut_type == "ins":
            if not (mutinfo[start_id] + 1 == mutinfo[end_id] and len(mutinfo[alt_id]) >= 1 and check_alt(
                    mutinfo[alt_id])):
                raise Exception("not in %s correct format!" % mut_type)
            if ref.fetch(mutinfo[chr_id], mutinfo[start_id] - 1, mutinfo[end_id]) == mutinfo[alt_id]:
                raise Exception("alt is same with reference!")
        elif mut_type == "del":
            if not ((mutinfo[alt_id] == "." or mutinfo[alt_id] == "") and mutinfo[end_id] >= mutinfo[start_id]):
                raise Exception("not in %s correct format!" % mut_type)
        elif mut_type == "sub":
            if not (len(mutinfo[alt_id]) >= 1 and mutinfo[end_id] >= mutinfo[start_id] and check_alt(mutinfo[alt_id])):
                raise Exception("not in %s correct format!" % mut_type)
            if ref.fetch(mutinfo[chr_id], mutinfo[start_id] - 1, mutinfo[end_id]) == mutinfo[alt_id]:
                raise Exception("alt is same with reference!")
        else:
            raise Exception("mutation type is not correct!")
        return True
    except Exception as e:
        return e.args


def check_polymorphism(bam, ref, mutinfo, maxsnpfrac):
    poly_statistic, reads_num = calPolymorphism(bam, mutinfo[0], mutinfo[1] - 1, mutinfo[2])
    # print poly_statistic, reads_num
    for pos in range(mutinfo[1] - 1, mutinfo[2]):
        ref_nul = ref.fetch(mutinfo[0], pos, pos + 1)
        ref_nul = ref_nul.upper()
        # print pos, ref_nul
        for nul in 'ATCG':
            if nul == ref_nul:
                continue
            if reads_num != 0:
                snp_frac = poly_statistic[pos][nul] * 1.0 / reads_num
            else:
                snp_frac = 0
            if snp_frac > maxsnpfrac:
                print snp_frac, maxsnpfrac
                return False
    return True


def check_mut_file(bam, ref, mut_file, snpfrac, invalid_log):
    chr_id, start_id, end_id, vaf_id, type_id, alt_id = (ind for ind in range(6))
    mut_list = []
    with open(mut_file, 'r') as fin_mut:
        for lineId, mutline in enumerate(fin_mut):
            if mutline[0] == "#":
                continue
            mutinfo = mutline.strip().split()
            # check mutation information
            flag = check_mut_info(mutinfo, ref)
            if not flag:
                invalid_log.info("line %s invalid in mutFile: %s" % (lineId, flag))
                continue
            # check single nucleotide polymorphism
            if mutinfo[type_id] in ["snv", "del", "sub"] and not check_polymorphism(bam, ref, mutinfo, snpfrac):
                invalid_log.info(
                    "haplotype in position %s:%s-%s: snp fraction bigger than run_args snvfrac(%s) " % (
                        mutinfo[chr_id], int(mutinfo[start_id])-1, int(mutinfo[end_id])-1, snpfrac))
                continue
            mut = Mutation(mutinfo[chr_id], mutinfo[start_id] - 1, mutinfo[end_id] - 1, mutinfo[vaf_id],
                           mutinfo[type_id], mutinfo[alt_id])
            mut_list.append(mut)
    fin_mut.close()
    return mut_list


def bed_cmp(mut_a, mut_b):
    if mut_a.chrom > mut_b.chrom:
        return 1
    elif mut_a.chrom == mut_b.chrom:
        if mut_a.start > mut_b.start:
            return 1
        elif mut_a.start == mut_b.start:
            if mut_a.end > mut_b.end:
                return 1
            else:
                return -1
        else:
            return -1
    else:
        return -1


def check_haplotype(bam, ref, mut_file, haplosize, snpfrac, invalid_log):
    mut_list = check_mut_file(bam, ref, mut_file, snpfrac, invalid_log)
    mut_list.sort(cmp=bed_cmp)
    cur_chrom, cur_end = None, None
    haplotype_list = []
    hap_mut_list = []
    for mut in mut_list:
        if not cur_end:
            cur_end = mut.end
            cur_chrom = mut.chrom
            hap_mut_list.append(mut)
            continue
        elif mut.chrom == cur_chrom and mut.start - cur_end < int(haplosize):
            if mut.start <= cur_end:
                invalid_log.info(
                    "haplotype in position %s:%s-%s: mutation invalid(overlapped by other mutations)" % (mut.chrom, mut.start, mut.end))
                continue
            hap_mut_list.append(mut)
        else:
            haplotype_list.append(hap_mut_list)
            hap_mut_list = list()
            hap_mut_list.append(mut)
        cur_chrom, cur_end = mut.chrom, mut.end
    haplotype_list.append(hap_mut_list)
    return haplotype_list


def get_haplotypes(bam_file, ref_fasta, mut_file, haplosize, snpfrac, invalid_log):
    ref = pysam.FastaFile(ref_fasta)
    bam = pysam.AlignmentFile(bam_file)
    haplotype_list_tmp = check_haplotype(bam, ref, mut_file, haplosize, snpfrac, invalid_log)
    haplotype_list = []
    for haplotypeMut in haplotype_list_tmp:
        haplotype_chrom, haplotype_start, haplotype_end, haplotype_freq = None, None, None, None
        for mut in haplotypeMut:
            if haplotype_chrom is None:
                haplotype_chrom = mut.chrom
            elif mut.chrom != haplotype_chrom:
                raise Exception("Error! This is a Bug")
            if haplotype_start is None or mut.start < haplotype_start:
                haplotype_start = mut.start
            if haplotype_end is None or mut.end > haplotype_end:
                haplotype_end = mut.end
            if haplotype_freq is None or mut.freq > haplotype_freq:
                haplotype_freq = mut.freq

        haplotype = Haplotype(haplotype_chrom, haplotype_start, haplotype_end, haplotype_freq, haplotypeMut)
        haplotype_list.append(haplotype)
    return haplotype_list
