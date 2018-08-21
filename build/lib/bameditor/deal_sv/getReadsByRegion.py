import os
from bameditor.common.bamconvert import bamIndex, getRegionReads


def get_reads_by_region(bam_file, sv_list, out_dir):
    # get reads by region bed
    region_bed_file = os.path.join(out_dir, "consider_region.bed")
    chrom_list = write_region_bed(region_bed_file, sv_list)
    exclude_bam_file_tmp = os.path.join(out_dir, "exclude_tmp.bam")
    used_bam_file_tmp = os.path.join(out_dir, "used_tmp.bam")
    getRegionReads(bam_file, region_bed_file, used_bam_file_tmp, exclude_bam_file_tmp)
    bamIndex(used_bam_file_tmp)
    return chrom_list, used_bam_file_tmp, exclude_bam_file_tmp


def write_region_bed(region_bed_file, sv_list):
    region_bed = open(region_bed_file, 'w')
    chrom_list = []
    for svinfo in sv_list:
        if svinfo.sv_type in ("dup", "cnv", "inv", "del"):
            if svinfo.chrom not in chrom_list:
                chrom_list.append(svinfo.chrom)
            region_bed.write("%s\t%s\t%s\n" % (svinfo.chrom, svinfo.start - 1000, svinfo.end + 1000))
        elif svinfo.sv_type in ("trans_balance", "trans_unbalance", "trans_chrom"):
            if svinfo.chrom not in chrom_list:
                chrom_list.append(svinfo.chrom)
            region_bed.write("%s\t%s\t%s\n" % (svinfo.chrom, svinfo.start - 1000, svinfo.end + 1000))
            if svinfo.trans_chrom not in chrom_list:
                chrom_list.append(svinfo.trans_chrom)
            region_bed.write(
                "%s\t%s\t%s\n" % (svinfo.trans_chrom, svinfo.trans_start - 1000, svinfo.trans_end + 1000))
    region_bed.close()
    return chrom_list
