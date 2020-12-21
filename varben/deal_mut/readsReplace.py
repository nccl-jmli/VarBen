from varben.common.methods import getReadStrand
from varben.common.bamconvert import remap, bamSort, bamIndex
import os
import pysam
from varben.common.methods import deal_life_reads, add_tag


def bam_add_tag(bam_file, out_bam_file):
    in_bam = pysam.AlignmentFile(bam_file, 'rb')
    out_bam = pysam.AlignmentFile(out_bam_file, 'wb', template=in_bam)
    for read in in_bam.fetch():
        read = add_tag(read)
        out_bam.write(read)
    in_bam.close()
    out_bam.close()


def reads_replace(bam_file, total_chosen_reads, seqer, flow_order, lib_key, barcode, tag, out_dir, aligner,
                  aligner_index, is_single):
    bam = pysam.AlignmentFile(bam_file)
    edit_bam_reads = {}
    for read in bam.fetch():
        read_name = read.query_name
        if read_name in total_chosen_reads:
            strand = getReadStrand(read)
            if read_name not in edit_bam_reads:
                edit_bam_reads[read_name] = {}
            if strand in total_chosen_reads[read_name]:
                my_read = total_chosen_reads[read_name][strand]
                read.query_sequence = my_read.query_sequence
                read.query_qualities = my_read.query_qualities
                if seqer == "life":
                    read = deal_life_reads(read, flow_order, lib_key, barcode)
                if tag:
                    read = add_tag(read)

                edit_bam_reads[read_name][strand] = read
            else:
                edit_bam_reads[read_name][strand] = read

    # write edited reads into edit.bam
    edit_bam_file = os.path.join(out_dir, "edit.bam")
    edit_bam = pysam.AlignmentFile(edit_bam_file, 'wb', template=bam)
    for read_name, readInfo in edit_bam_reads.items():
        for strand, read in readInfo.items():
            edit_bam.write(read)
    edit_bam.close()

    # write not edited reads into exclude.bam
    exclude_bam_file = os.path.join(out_dir, "exclude.bam")
    exclude_bam = pysam.AlignmentFile(exclude_bam_file, 'wb', template=bam)
    for read in bam.fetch():
        read_name = read.query_name
        if read_name not in edit_bam_reads:
            exclude_bam.write(read)
    exclude_bam.close()

    # remap the edited reads
    header = os.path.join(out_dir, 'bam.header')
    os.system('samtools view -H %s|grep "^@RG" > %s' % (bam_file, header))
    head = open(header, 'r').readline().rstrip().replace('\t','\\t')
    if not head:
        head = None
    edit_remap_bam_file = os.path.join(out_dir, "edit.remap.bam")
    remap(aligner_index, edit_bam_file, edit_remap_bam_file, aligner, is_single, header=head)
    edit_remap_bam_sorted_prefix = os.path.join(out_dir, "edit.remap.sort")
    edit_remap_bam_sorted_file = os.path.join(out_dir, "edit.remap.sort.bam")
    bamSort(edit_remap_bam_file, edit_remap_bam_sorted_prefix)
    bamIndex(edit_remap_bam_sorted_file)
    if tag:
        edit_remap_addtag_file = os.path.join(out_dir, "edit.remap.sort.bam")
        bam_add_tag(edit_remap_bam_sorted_file, edit_remap_addtag_file)
    else:
        edit_remap_addtag_file = edit_remap_bam_sorted_file

    return edit_remap_addtag_file, exclude_bam_file


