import os
import pysam
import copy
from bameditor.common.methods import getKeyName, get_new_readname
from bameditor.common.methods import deal_life_reads, add_tag
from bameditor.common.bamconvert import remap, bamAddRG, bamIndex
from array import array


def get_newname_dict(readname_list_file):
    modify_read_name_dict = {}
    with open(readname_list_file, 'r') as readname_list:
        for line in readname_list:
            reads = line.strip().split(",")
            new_name = get_new_readname()
            for readname in reads:
                modify_read_name_dict[readname] = new_name
    readname_list.close()
    return modify_read_name_dict


def get_sequence_dict(total_reads_file):
    seq_dict = {}
    quan_dict = {}
    with open(total_reads_file, 'r') as fin:
        for line in fin:
            if not line:
                break
            data = line.strip().split("\t")
            seq_dict[data[0]] = data[2]
            s = data[3].index('[')
            e = data[3].index(']')
            quan = eval(data[3][s:e+1])
            # print s, e, quan
            quan_dict[data[0]] = array('B', quan)
    fin.close()
    return seq_dict, quan_dict


def merge_edit_bam(bam_file, out_dir, is_single, total_modify_reads_file, total_add_reads_file, used_bam_file,
                   total_modify_list, total_add_list, seqer, aligner,
                   aligner_index, flow_order, lib_key, barcode, tag):
    bam = pysam.AlignmentFile(bam_file, 'rb')
    edit_bam_file = os.path.join(out_dir, "edit.bam")
    edit_bam = pysam.AlignmentFile(edit_bam_file, 'wb', template=bam)
    readname_convert_file = os.path.join(out_dir, "readname_convert.txt")
    fout_convert = open(readname_convert_file, 'w')
    # edit_bam_reads = {}
    used_bam = pysam.AlignmentFile(used_bam_file, 'rb')
    used_reads = {}
    for read in used_bam.fetch():
        keyname = getKeyName(read)
        used_reads[keyname] = read
    used_bam.close()
    # modify_read_name_dict = get_newname_dict(total_modify_list)
    # add_read_name_dict = get_newname_dict(total_add_list)

    modify_reads_seq, modify_reads_quan = get_sequence_dict(total_modify_reads_file)
    add_reads_seq, add_reads_quan = get_sequence_dict(total_add_reads_file)
    if is_single:
        with open(total_modify_list) as fin:
            for line in fin:
                if not line:
                    break
                data = line.strip().split(",")
                read1_name = data[0]
                if read1_name not in used_reads:
                    continue
                orig_read1 = used_reads[read1_name]
                new_read1 = copy.deepcopy(orig_read1)
                new_read1.query_sequence = modify_reads_seq[read1_name]
                new_read1.query_qualities = modify_reads_quan[read1_name]
                new_name = get_new_readname()
                new_read1.query_name = new_name
                if seqer == "life":
                    new_read1 = deal_life_reads(new_read1, flow_order, lib_key, barcode)
                fout_convert.write("%s: %s, %s, %s-%s\n" % (
                    new_name, orig_read1.query_name, new_read1.is_read1, new_read1.reference_start,
                    new_read1.reference_end))
                edit_bam.write(new_read1)
        fin.close()
        with open(total_add_list) as fin:
            for line in fin:
                if not line:
                    break
                data = line.strip().split(",")
                read1_name = data[0]
                if read1_name not in used_reads:
                    continue
                orig_read1 = used_reads[read1_name]
                new_read1 = copy.deepcopy(orig_read1)
                new_read1.query_sequence = add_reads_seq[read1_name]
                new_read1.query_qualities = add_reads_quan[read1_name]
                new_name = get_new_readname()
                new_read1.query_name = new_name
                fout_convert.write("%s: %s, %s, %s-%s\n" % (
                    new_name, orig_read1.query_name, new_read1.is_read1, new_read1.reference_start,
                    new_read1.reference_end))
                edit_bam.write(new_read1)
        fin.close()

    else:
        with open(total_modify_list) as fin:
            for line in fin:
                if not line:
                    break
                data = line.strip().split(",")
                read1_name, read2_name = data[0], data[1]
                if read1_name not in used_reads or read2_name not in used_reads:
                    continue
                orig_read1 = used_reads[read1_name]
                orig_read2 = used_reads[read2_name]
                new_read1 = copy.deepcopy(orig_read1)
                new_read2 = copy.deepcopy(orig_read2)
                print read1_name
                new_read1.query_sequence = modify_reads_seq[read1_name]
                new_read1.query_qualities = modify_reads_quan[read1_name]
                new_read2.query_sequence = modify_reads_seq[read2_name]
                new_read2.query_qualities = modify_reads_quan[read2_name]
                new_name = get_new_readname()
                new_read1.query_name = new_name
                new_read2.query_name = new_name
                fout_convert.write("%s: %s, %s, %s-%s\n" % (
                            new_name, orig_read1.query_name, new_read1.is_read1, new_read1.reference_start,
                            new_read1.reference_end))
                fout_convert.write("%s: %s, %s, %s-%s\n" % (
                    new_name, orig_read2.query_name, new_read2.is_read1, new_read2.reference_start,
                    new_read2.reference_end))
                edit_bam.write(new_read1)
                edit_bam.write(new_read2)
        fin.close()
        with open(total_add_list) as fin:
            for line in fin:
                if not line:
                    break
                data = line.strip().split(",")
                read1_name, read2_name = data[0], data[1]
                if read1_name not in used_reads or read2_name not in used_reads:
                    continue
                orig_read1 = used_reads[read1_name]
                orig_read2 = used_reads[read2_name]
                new_read1 = copy.deepcopy(orig_read1)
                new_read2 = copy.deepcopy(orig_read2)
                new_read1.query_sequence = add_reads_seq[read1_name]
                new_read1.query_qualities = add_reads_quan[read1_name]
                new_read2.query_sequence = add_reads_seq[read2_name]
                new_read2.query_qualities = add_reads_quan[read2_name]
                new_name = get_new_readname()
                new_read1.query_name = new_name
                new_read2.query_name = new_name
                fout_convert.write("%s: %s, %s, %s-%s\n" % (
                            new_name, orig_read1.query_name, new_read1.is_read1, new_read1.reference_start,
                            new_read1.reference_end))
                fout_convert.write("%s: %s, %s, %s-%s\n" % (
                    new_name, orig_read2.query_name, new_read2.is_read1, new_read2.reference_start,
                    new_read2.reference_end))
                edit_bam.write(new_read1)
                edit_bam.write(new_read2)
        fin.close()
    edit_bam.close()
    edit_remap_bam_file = os.path.join(out_dir, "edit.remap.bam")
    remap(aligner_index, edit_bam_file, edit_remap_bam_file, aligner, is_single)
    return edit_remap_bam_file

    # for read in modify_reads.fetch():
    #     old_name = read.query_name
    #     read_name = getKeyName(read)
    #     new_read1 = copy.deepcopy(read)
    #     new_read1.query_name = modify_read_name_dict[read_name]
    #     edit_bam.write(new_read1)
    # for read in add_reads.fetch():
    #     read_name = getKeyName(read)
    #     new_read1 = copy.deepcopy(read)
    #     new_read1.query_name = add_read_name_dict[read_name]
    #     edit_bam.write(new_read1)
    # edit_bam.close()


    # if is_single:
    #     pass
        # for read_pair in total_modify_reads:
        #     read1 = read_pair[0]
        #     keyname_read1 = getKeyName(read1)
        #     orig_read1 = used_reads[keyname_read1]
        #     new_read1 = copy.deepcopy(orig_read1)
        #     new_read1.query_sequence = read1.query_sequence
        #     new_read1.query_qualities = read1.query_qualities
        #     new_name = read1.query_name.split(":")[0] + ":" + get_new_readname()
        #     new_read1.query_name = new_name
        #     if seqer == "life":
        #         new_read1 = deal_life_reads(new_read1, flow_order, lib_key, barcode)
        #     if tag:
        #         new_read1 = add_tag(new_read1)
        #     edit_bam.write(new_read1)
        #     fout_convert.write("%s: %s, %s, %s-%s\n" % (
        #         new_name, orig_read1.query_name, new_read1.is_read1, new_read1.reference_start,
        #         new_read1.reference_end))
        #     strand = getReadStrand(new_read1)
        #     if new_name not in edit_bam_reads:
        #         edit_bam_reads[new_name] = dict()
        #     edit_bam_reads[new_name][strand] = new_read1
        #
        # for read_pair in total_add_reads:
        #     read1 = read_pair[0]
        #     keyname_read1 = getKeyName(read1)
        #     orig_read1 = used_reads[keyname_read1]
        #     new_read1 = copy.deepcopy(orig_read1)
        #     new_name = get_new_readname()
        #     new_read1.query_name = new_name
        #     edit_bam.write(new_read1)
        #     fout_convert.write("%s: %s, %s, %s-%s\n" % (
        #         new_name, orig_read1.query_name, new_read1.is_read1, new_read1.reference_start,
        #         new_read1.reference_end))
        #     strand = getReadStrand(new_read1)
        #     if new_name not in edit_bam_reads:
        #         edit_bam_reads[new_name] = dict()
        #     edit_bam_reads[new_name][strand] = new_read1

    # else:
    #     for read_pair in total_modify_reads + total_add_reads:
    #         read1 = read_pair[0]
    #         read2 = read_pair[1]
    #         keyname_read1 = getKeyName(read1)
    #         keyname_read2 = getKeyName(read2)
    #         orig_read1 = used_reads[keyname_read1]
    #         orig_read2 = used_reads[keyname_read2]
    #         orig_read1_name = orig_read1.query_name
    #         orig_read2_name = orig_read2.query_name
    #         new_read1 = copy.deepcopy(orig_read1)
    #         new_read2 = copy.deepcopy(orig_read2)
    #         new_read1.query_sequence = read1.query_sequence
    #         new_read1.query_qualities = read1.query_qualities
    #         new_read2.query_sequence = read2.query_sequence
    #         new_read2.query_qualities = read2.query_qualities
    #         new_name = get_new_readname()
    #         new_read1.query_name = new_name
    #         new_read2.query_name = new_name
    #         strand1 = getReadStrand(new_read1)
    #         strand2 = getReadStrand(new_read2)
    #         if new_name not in edit_bam_reads:
    #             edit_bam_reads[new_name] = dict()
    #         edit_bam_reads[new_name][strand1] = new_read1
    #         edit_bam_reads[new_name][strand2] = new_read2
    #         if tag:
    #             new_read1 = add_tag(new_read1)
    #             new_read2 = add_tag(new_read2)
    #
    #         fout_convert.write("%s: %s, %s, %s, %s, %s-%s, %s-%s\n" % (
    #             new_name, orig_read1_name, orig_read2_name, new_read1.is_read1, new_read2.is_read2,
    #             new_read1.reference_start, new_read1.reference_end, new_read2.reference_start,
    #             new_read2.reference_end,))
    #         edit_bam.write(new_read1)
    #         edit_bam.write(new_read2)
    # fout_convert.close()
    # edit_bam.close()

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
