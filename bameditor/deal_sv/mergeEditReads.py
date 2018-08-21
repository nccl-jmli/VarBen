from bameditor.common.methods import getKeyName
import os


def _get_write(total_reads, reads_file_out, reads_pair_out):
    reads_file = open(reads_file_out, 'w')
    reads_pair = open(reads_pair_out, 'w')
    for read_pair in total_reads:
        tmp = []
        for read in read_pair:
            read_name = getKeyName(read)
            tmp.append(read_name)
            reads_file.write(str(read))
        reads_pair.write("%s\n" % ",".join(tmp))
    reads_file.close()
    reads_pair.close()


def get_write_reads(total_modify_reads, total_add_reads, total_delete_reads, bamfile, out_dir):
    reads_file_list = []
    reads_pair_list = []
    for typ, reads_dict in zip(('modify', 'add', 'delete'), (total_modify_reads, total_add_reads, total_delete_reads)):
        reads_file = os.path.join(out_dir, "total_%s_reads.txt" % typ)
        pair_list = os.path.join(out_dir, "total_%s_pair.list" % typ)

        _get_write(reads_dict, reads_file, pair_list)
        reads_pair_list.append(pair_list)
        reads_file_list.append(reads_file)
    return reads_file_list + reads_pair_list

