from bameditor.common.methods import getKeyName


def merge_edit_reads(total_modify_reads, total_add_reads, total_delete_reads):
    total_modify_readname_list, total_add_readname_list, total_delete_readname_list = [], [], []
    total_modify_reads_list, total_add_reads_list, total_delete_reads_list = [], [], []

    for read_pair in total_modify_reads:
        for read in read_pair:
            read_name = getKeyName(read)
            total_modify_readname_list.append(read_name)
            total_modify_reads_list.append(read)
    for read_pair in total_add_reads:
        for read in read_pair:
            read_name = getKeyName(read)
            total_add_readname_list.append(read_name)
            total_add_reads_list.append(read)
    for read_pair in total_delete_reads:
        for read in read_pair:
            read_name = getKeyName(read)
            total_delete_readname_list.append(read_name)
            total_delete_reads_list.append(read)

    return total_modify_readname_list, total_delete_readname_list, total_add_readname_list
