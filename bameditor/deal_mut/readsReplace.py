from bameditor.common.methods import getReadStrand
from bameditor.common.bamconvert import remap, bamAddRG
import os
import pysam
from bameditor.common.methods import deal_life_reads, add_tag


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
    edit_remap_bam_file = os.path.join(out_dir, "edit.remap.bam")
    remap(aligner_index, edit_bam_file, edit_remap_bam_file, aligner, is_single)

    # if not is_single:
    #     edit_remap = pysam.AlignmentFile(edit_remap_bam_file, 'rb')
    #     editRemapBam_addRG_File = os.path.join(out_dir, "edit.remap.addRG.bam")
    #     bamAddRG(edit_remap, edit_bam_reads, bam, editRemapBam_addRG_File)
    #     edit_remap.close()
    # else:
    #     editRemapBam_addRG_File = edit_remap_bam_file

    return edit_remap_bam_file, exclude_bam_file


