from bameditor.deal_mut.readEditor import editRead
from multiprocessing import Pool


def reads_modify(total_chosen_reads, total_chosen_reads_muts, ref_fasta, prceoss):
    edit_pool = Pool(processes=int(prceoss))
    edit_res = []
    keys_list = []
    for read_name in total_chosen_reads:
        for strand in total_chosen_reads[read_name]:
            keys_list.append((read_name, strand))
            read = total_chosen_reads[read_name][strand]
            mut_list = total_chosen_reads_muts[read_name][strand]
            edit_res.append(edit_pool.apply_async(editRead, args=(read, ref_fasta, mut_list)))
    edit_pool.close()
    edit_pool.join()

    # step5: write edited reads to fastq file
    for keys, res_tmp in zip(keys_list, edit_res):
        res = res_tmp.get()
        if not res:
            continue
        sequence, quality, shift = res
        read_name, strand = keys
        read = total_chosen_reads[read_name][strand]
        read.query_sequence = sequence
        read.query_qualities = quality
    return True
