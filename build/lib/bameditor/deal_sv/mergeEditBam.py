import os
import pysam
import copy
from bameditor.common.methods import getKeyName, get_new_readname
from bameditor.common.methods import deal_life_reads, add_tag
from bameditor.common.bamconvert import remap, bamAddRG, bamIndex
from bameditor.common.methods import getReadStrand


def merge_edit_bam(bam_file, out_dir, is_single, total_modify_reads, total_add_reads, used_reads, seqer, aligner,
                   aligner_index, flow_order, lib_key, barcode, tag):
    bam = pysam.AlignmentFile(bam_file, 'rb')
    edit_bam_file = os.path.join(out_dir, "edit.bam")
    edit_bam = pysam.AlignmentFile(edit_bam_file, 'wb', template=bam)
    readname_convert_file = os.path.join(out_dir, "readname_convert.txt")
    fout_convert = open(readname_convert_file, 'w')
    edit_bam_reads = {}
    if is_single:
        for read_pair in total_modify_reads:
            read1 = read_pair[0]
            keyname_read1 = getKeyName(read1)
            orig_read1 = used_reads[keyname_read1]
            new_read1 = copy.deepcopy(orig_read1)
            new_read1.query_sequence = read1.query_sequence
            new_read1.query_qualities = read1.query_qualities
            new_name = read1.query_name.split(":")[0] + ":" + get_new_readname()
            new_read1.query_name = new_name
            if seqer == "life":
                new_read1 = deal_life_reads(new_read1, flow_order, lib_key, barcode)
            if tag:
                new_read1 = add_tag(new_read1)
            edit_bam.write(new_read1)
            fout_convert.write("%s: %s, %s, %s-%s\n" % (
                new_name, orig_read1.query_name, new_read1.is_read1, new_read1.reference_start,
                new_read1.reference_end))
            strand = getReadStrand(new_read1)
            if new_name not in edit_bam_reads:
                edit_bam_reads[new_name] = dict()
            edit_bam_reads[new_name][strand] = new_read1

        for read_pair in total_add_reads:
            read1 = read_pair[0]
            keyname_read1 = getKeyName(read1)
            orig_read1 = used_reads[keyname_read1]
            new_read1 = copy.deepcopy(orig_read1)
            new_name = get_new_readname()
            new_read1.query_name = new_name
            edit_bam.write(new_read1)
            fout_convert.write("%s: %s, %s, %s-%s\n" % (
                new_name, orig_read1.query_name, new_read1.is_read1, new_read1.reference_start,
                new_read1.reference_end))
            strand = getReadStrand(new_read1)
            if new_name not in edit_bam_reads:
                edit_bam_reads[new_name] = dict()
            edit_bam_reads[new_name][strand] = new_read1

    else:
        for read_pair in total_modify_reads + total_add_reads:
            read1 = read_pair[0]
            read2 = read_pair[1]
            keyname_read1 = getKeyName(read1)
            keyname_read2 = getKeyName(read2)
            orig_read1 = used_reads[keyname_read1]
            orig_read2 = used_reads[keyname_read2]
            orig_read1_name = orig_read1.query_name
            orig_read2_name = orig_read2.query_name
            new_read1 = copy.deepcopy(orig_read1)
            new_read2 = copy.deepcopy(orig_read2)
            new_read1.query_sequence = read1.query_sequence
            new_read1.query_qualities = read1.query_qualities
            new_read2.query_sequence = read2.query_sequence
            new_read2.query_qualities = read2.query_qualities
            new_name = get_new_readname()
            new_read1.query_name = new_name
            new_read2.query_name = new_name
            strand1 = getReadStrand(new_read1)
            strand2 = getReadStrand(new_read2)
            if new_name not in edit_bam_reads:
                edit_bam_reads[new_name] = dict()
            edit_bam_reads[new_name][strand1] = new_read1
            edit_bam_reads[new_name][strand2] = new_read2
            if tag:
                new_read1 = add_tag(new_read1)
                new_read2 = add_tag(new_read2)

            fout_convert.write("%s: %s, %s, %s, %s, %s-%s, %s-%s\n" % (
                new_name, orig_read1_name, orig_read2_name, new_read1.is_read1, new_read2.is_read2,
                new_read1.reference_start, new_read1.reference_end, new_read2.reference_start,
                new_read2.reference_end,))
            edit_bam.write(new_read1)
            edit_bam.write(new_read2)
    fout_convert.close()
    edit_bam.close()

    edit_remap_bam_file = os.path.join(out_dir, "edit.remap.bam")
    remap(aligner_index, edit_bam_file, edit_remap_bam_file, aligner, is_single)

    if not is_single:
        editRemap = pysam.AlignmentFile(edit_remap_bam_file, 'rb')
        editRemapBam_addRG_File = os.path.join(out_dir, "edit.remap.addRG.bam")
        bamAddRG(editRemap, edit_bam_reads, bam, editRemapBam_addRG_File)
        editRemap.close()
    else:
        editRemapBam_addRG_File = edit_remap_bam_file
    bamIndex(editRemapBam_addRG_File)
    return editRemapBam_addRG_File
