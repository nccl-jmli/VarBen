import os
import pysam
import copy
from varben.common.methods import getKeyName, get_new_readname
from varben.common.methods import deal_life_reads
from varben.common.bamconvert import remap, bamAddRG, bamIndex
from array import array
from varben.common.bamconvert import bamSort


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

    header = os.path.join(out_dir, 'bam.header')
    os.system('samtools view -H %s|grep "^@RG" > %s' % (bam_file, header))
    head = open(header, 'r').readline().rstrip().replace('\t','\\t')
    if not head:
        head = None
    edit_remap_bam_file = os.path.join(out_dir, "edit.remap.bam")
    remap(aligner_index, edit_bam_file, edit_remap_bam_file, aligner, is_single, head)
    edit_remap_bam_sorted_prefix = os.path.join(out_dir, "edit.remap.sort")
    edit_remap_bam_sorted_file = os.path.join(out_dir, "edit.remap.sort.bam")
    bamSort(edit_remap_bam_file, edit_remap_bam_sorted_prefix)
    return edit_remap_bam_sorted_file
