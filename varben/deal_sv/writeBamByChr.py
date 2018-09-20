__author__ = 'fangshuangsang'
import pysam
from varben.common.methods import getKeyName
import os
from multiprocessing.pool import Pool
from varben.common.bamconvert import bamIndex, bamMerge, bamSort


def get_name_list(reads_file, chr):
    fin_reads = open(reads_file, 'r')
    total_readname_list = []
    for read_info in fin_reads:
        if not read_info:
            break
        tmp = read_info.strip().split("\t")
        # print tmp, chr
        if tmp[1] == chr:
            total_readname_list.append(tmp[0])
    return total_readname_list


def write_bam_byChr(bamFile, chr, excludeBamFile, editBamFile, total_modify_reads_file, total_delete_reads_file,
                    total_add_reads_file):
    print chr
    bam = pysam.AlignmentFile(bamFile, 'rb')
    excludeBam = pysam.AlignmentFile(excludeBamFile, 'wb', template=bam)
    editBam = pysam.AlignmentFile(editBamFile, 'wb', template=bam)
    delete = open(editBamFile + ".del", 'w')
    m = 0
    modifyReadsName = get_name_list(total_modify_reads_file, chr)
    deleteReadsName = get_name_list(total_delete_reads_file, chr)
    addReadsName = get_name_list(total_add_reads_file, chr)
    for read in bam.fetch(chr):
        m += 1
        keyname = getKeyName(read)
        if keyname in modifyReadsName:
            editBam.write(read)
        elif keyname in deleteReadsName:
            delete.write(keyname + "\n")
            continue
        elif keyname in addReadsName:
            editBam.write(read)
            excludeBam.write(read)
        else:
            excludeBam.write(read)
    delete.close()
    print "Total reads: ", m
    bam.close()
    excludeBam.close()
    editBam.close()


def write_sub_bam(chrom_list, used_bam_file_tmp, exclude_bam_file_tmp, out_dir, total_modify_reads_file,
                  total_delete_reads_file, total_add_reads_file, process):
    write_bam_pool = Pool(int(process))
    exclude_bam_list = [exclude_bam_file_tmp]
    usedBamList = []
    for chrom in chrom_list:
        excludeBam_chr = "%s/exclude_%s.bam" % (out_dir, chrom)
        exclude_bam_list.append(excludeBam_chr)
        usedBam_chr = "%s/used_%s.bam" % (out_dir, chrom)
        usedBamList.append(usedBam_chr)

        write_bam_pool.apply_async(write_bam_byChr, args=(used_bam_file_tmp, chrom, excludeBam_chr, usedBam_chr,
                                                          total_modify_reads_file, total_delete_reads_file,
                                                          total_add_reads_file))
    write_bam_pool.close()
    write_bam_pool.join()

    exclude_bam_file = os.path.join(out_dir, "exclude.bam")
    bamMerge(exclude_bam_list, exclude_bam_file)
    used_bam_file = os.path.join(out_dir, "used.bam")
    if len(usedBamList) != 1:
        bamMerge(usedBamList, used_bam_file)
    else:
        used_bam_file = usedBamList[0]

    bamSort(used_bam_file, os.path.join(out_dir, "used.sort"))
    used_sort_bam_file = os.path.join(out_dir, "used.sort.bam")
    bamIndex(used_sort_bam_file)
    return used_sort_bam_file, exclude_bam_file
