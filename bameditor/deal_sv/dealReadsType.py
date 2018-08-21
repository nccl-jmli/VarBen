import random
import re
import copy
from bameditor.common.methods import check_reads_pair, getComplementarySeq

try_max_time = 100
shift_num = 1000


def mend_read_part(ref, read, start, end, readsType, relPart, svtype="del",
                   subPos=None):  # relPart: which has existed, for example: "left" means left existed, need mend right part
    seq_len = len(read.query_sequence)
    posPairList = read.get_aligned_pairs
    svPosIndex = 0
    if readsType == "type1" or readsType == "type6":
        if relPart == "left":
            for i, pair in enumerate(posPairList):
                svPosIndex = pair[0]
                if pair[1] >= start:
                    break
            if svtype == "del":
                pos_start = end + 1
                pos_end = pos_start + seq_len - svPosIndex
                mend_seq = ref.fetch(read.reference_name, pos_start, pos_end)
                new_seq = read.query_sequence[:svPosIndex] + mend_seq
            elif svtype == "inv":
                pos_end = end + 1
                pos_start = pos_end - (seq_len - svPosIndex)
                mend_seq = ref.fetch(read.reference_name, pos_start, pos_end)
                new_seq = read.query_sequence[:svPosIndex] + getComplementarySeq(mend_seq)[::-1]
            elif svtype == "trans_balance" or svtype == "trans_chrom":
                res = re.match("(\w*):(\d*)-(\d*)", subPos)
                trans_chr, trans_start, trans_end = res.group(1), int(res.group(2)), int(res.group(3))
                pos_start = trans_start
                pos_end = pos_start + (seq_len - svPosIndex)
                mend_seq = ref.fetch(trans_chr, pos_start, pos_end)
                new_seq = read.query_sequence[:svPosIndex] + mend_seq
            elif svtype == "trans_unbalance":
                res = re.match("(\w*):(\d*)-(\d*)", subPos)
                trans_chr, trans_start, trans_end = res.group(1), int(res.group(2)), int(res.group(3))
                pos_start = trans_start
                pos_end = pos_start + seq_len - (svPosIndex + 1)
                mend_seq = ref.fetch(trans_chr, pos_start, pos_end)
                new_seq = read.query_sequence[:svPosIndex + 1] + mend_seq
            elif readsType == "type6" and svtype == "dup":
                pos_end = end + 1
                pos_start = pos_end - svPosIndex
                mend_seq = ref.fetch(read.reference_name, pos_start, pos_end)
                new_seq = mend_seq + read.query_sequence[svPosIndex:]

        else:
            for i, pair in enumerate(posPairList):
                svPosIndex = pair[0]
                if pair[1] >= end:
                    break
            if svtype == "del":
                pos_end = start - 1 + 1
                pos_start = pos_end - (svPosIndex + 1)
                mend_seq = ref.fetch(read.reference_name, pos_start, pos_end)
                new_seq = mend_seq + read.query_sequence[svPosIndex + 1:]
            elif svtype == "inv":
                pos_start = start
                pos_end = start + (svPosIndex + 1)
                mend_seq = ref.fetch(read.reference_name, pos_start, pos_end)
                new_seq = getComplementarySeq(mend_seq)[::-1] + read.query_sequence[svPosIndex + 1:]
            elif svtype == "trans_balance":
                res = re.match("(\w*):(\d*)-(\d*)", subPos)
                trans_chr, trans_start, trans_end = res.group(1), int(res.group(2)), int(res.group(3))
                pos_end = trans_end + 1
                pos_start = pos_end - (svPosIndex + 1)
                mend_seq = ref.fetch(trans_chr, pos_start, pos_end)
                new_seq = mend_seq + read.query_sequence[svPosIndex + 1:]
            elif readsType == "type6" and svtype == "dup":
                pos_start = start
                pos_end = pos_start + seq_len - (svPosIndex + 1)
                mend_seq = ref.fetch(read.reference_name, pos_start, pos_end)
                new_seq = read.query_sequence[:svPosIndex + 1] + mend_seq
            elif readsType == "type6" and svtype == "trans_unbalance":
                res = re.match("(\w*):(\d*)-(\d*)", subPos)
                trans_chr, trans_start, trans_end = res.group(1), int(res.group(2)), int(res.group(3))
                pos_end = end + 1
                pos_start = pos_end - (svPosIndex + 1)
                mend_seq = ref.fetch(trans_chr, pos_start, pos_end)
                new_seq = mend_seq + read.query_sequence[svPosIndex + 1:]

    elif readsType == "type2":
        if relPart == "left":
            for i, pair in enumerate(posPairList):
                svPosIndex = pair[0]
                if pair[1] >= start:
                    break
            if svtype == "inv":
                pos_start = end + 1
                pos_end = pos_start + svPosIndex
                mend_seq = ref.fetch(read.reference_name, pos_start, pos_end)
                new_seq = getComplementarySeq(mend_seq)[::-1] + read.query_sequence[svPosIndex:]
            elif svtype == "dup":
                pos_end = end + 1
                pos_start = pos_end - svPosIndex
                mend_seq = ref.fetch(read.reference_name, pos_start, pos_end)
                new_seq = mend_seq + read.query_sequence[svPosIndex:]

            elif svtype == "trans_balance" or svtype == "trans_chrom":
                res = re.match("(\w*):(\d*)-(\d*)", subPos)
                trans_chr, trans_start, trans_end = res.group(1), int(res.group(2)), int(res.group(3))
                pos_end = trans_start
                pos_start = pos_end - svPosIndex
                mend_seq = ref.fetch(trans_chr, pos_start, pos_end)
                new_seq = mend_seq + read.query_sequence[svPosIndex:]

        elif relPart == "right":
            for i, pair in enumerate(posPairList):
                svPosIndex = pair[0]
                if pair[1] >= end:
                    break
            if svtype == "inv":
                pos_end = start
                pos_start = pos_end - (seq_len - svPosIndex - 1)
                mend_seq = ref.fetch(read.reference_name, pos_start, pos_end)

                new_seq = read.query_sequence[:svPosIndex + 1] + getComplementarySeq(mend_seq)[::-1]
                # print read.query_name, read.query_sequence, mend_seq, new_seq
            elif svtype == "dup":
                pos_start = start
                pos_end = pos_start + seq_len - (svPosIndex + 1)
                mend_seq = ref.fetch(read.reference_name, pos_start, pos_end)
                new_seq = read.query_sequence[0:svPosIndex + 1] + mend_seq
            elif svtype == "trans_balance":
                res = re.match("(\w*):(\d*)-(\d*)", subPos)
                trans_chr, trans_start, trans_end = res.group(1), int(res.group(2)), int(res.group(3))
                pos_start = trans_end + 1
                pos_end = pos_start + seq_len - (svPosIndex + 1)
                mend_seq = ref.fetch(trans_chr, pos_start, pos_end)
                new_seq = read.query_sequence[:svPosIndex + 1] + mend_seq
            elif svtype == "trans_unbalance":
                res = re.match("(\w*):(\d*)-(\d*)", subPos)
                trans_chr, trans_start, trans_end = res.group(1), int(res.group(2)), int(res.group(3))
                pos_end = trans_end + 1
                pos_start = pos_end - svPosIndex
                mend_seq = ref.fetch(trans_chr, pos_start, pos_end)
                # print trans_chr, pos_start, pos_end
                new_seq = mend_seq + read.query_sequence[svPosIndex:]

    elif readsType == "type5":
        if relPart == "left":
            for i, pair in enumerate(posPairList):
                svPosIndex = pair[0]
                if pair[1] >= start:
                    break
            if svtype == "del":
                pos_start = end + 1
                pos_end = pos_start + seq_len - svPosIndex
                mend_seq = ref.fetch(read.reference_name, pos_start, pos_end)
                new_seq = read.query_sequence[:svPosIndex] + mend_seq

            elif svtype == "dup":
                pos_end = end + 1
                pos_start = pos_end - svPosIndex
                mend_seq = ref.fetch(read.reference_name, pos_start, pos_end)
                new_seq = mend_seq + read.query_sequence[svPosIndex:]

            elif svtype == "trans_balance" or svtype == "trans_chrom":
                res = re.match("(\w*):(\d*)-(\d*)", subPos)
                trans_chr, trans_start, trans_end = res.group(1), int(res.group(2)), int(res.group(3))
                pos_start = trans_start
                pos_end = pos_start + (seq_len - svPosIndex)
                mend_seq = ref.fetch(trans_chr, pos_start, pos_end)
                new_seq = read.query_sequence[:svPosIndex] + mend_seq
            elif svtype == "trans_unbalance":
                res = re.match("(\w*):(\d*)-(\d*)", subPos)
                trans_chr, trans_start, trans_end = res.group(1), int(res.group(2)), int(res.group(3))
                pos_start = trans_start
                pos_end = pos_start + seq_len - (svPosIndex + 1)
                mend_seq = ref.fetch(trans_chr, pos_start, pos_end)
                new_seq = read.query_sequence[:svPosIndex + 1] + mend_seq

            elif svtype == "inv":
                pos_end = end + 1
                pos_start = pos_end - (seq_len - svPosIndex)
                mend_seq = ref.fetch(read.reference_name, pos_start, pos_end)
                new_seq = read.query_sequence[:svPosIndex] + getComplementarySeq(mend_seq)[::-1]

        elif relPart == "right":
            for i, pair in enumerate(posPairList):
                svPosIndex = pair[0]
                if pair[1] >= end:
                    break

            if svtype == "dup":
                pos_start = start
                pos_end = pos_start + seq_len - (svPosIndex + 1)
                mend_seq = ref.fetch(read.reference_name, pos_start, pos_end)
                new_seq = read.query_sequence[0:svPosIndex + 1] + mend_seq
            elif svtype == "trans_balance":
                res = re.match("(\w*):(\d*)-(\d*)", subPos)
                trans_chr, trans_start, trans_end = res.group(1), int(res.group(2)), int(res.group(3))
                pos_end = trans_end + 1
                pos_start = pos_end - (svPosIndex + 1)
                mend_seq = ref.fetch(trans_chr, pos_start, pos_end)
                new_seq = mend_seq + read.query_sequence[svPosIndex + 1:]
            elif svtype == "trans_unbalance":
                res = re.match("(\w*):(\d*)-(\d*)", subPos)
                trans_chr, trans_start, trans_end = res.group(1), int(res.group(2)), int(res.group(3))
                pos_end = trans_end + 1
                pos_start = pos_end - svPosIndex
                mend_seq = ref.fetch(trans_chr, pos_start, pos_end)
                new_seq = mend_seq + read.query_sequence[svPosIndex:]
            elif svtype == "inv":
                pos_start = start
                pos_end = start + (svPosIndex + 1)
                mend_seq = ref.fetch(read.reference_name, pos_start, pos_end)
                new_seq = getComplementarySeq(mend_seq)[::-1] + read.query_sequence[svPosIndex + 1:]

    # if len(new_seq) != seq_len:
    #     print readsType, relPart, svtype
    return new_seq


def deal_type1(ref, reads_type1_left, reads_type1_right, freq, start, end, svtype, subPos=None):
    # fix up
    print "deal type1 start......"
    reads_left_num = len(reads_type1_left)
    reads_right_num = len(reads_type1_right)
    reads_left_mend_id = random_mendIDList(reads_left_num, freq)
    reads_right_mend_id = random_mendIDList(reads_right_num, freq)
    total_reads = []
    for read_pair_id in reads_left_mend_id:
        read_pair = reads_type1_left[read_pair_id]
        read = read_pair[1]
        read_mate = read_pair[0]
        new_seq = mend_read_part(ref, read, start, end, "type1", "left", svtype, subPos)
        qual = read.query_qualities
        read.query_sequence = new_seq
        read.query_qualities = qual
        total_reads.append([read_mate, read])

    for read_pair_id in reads_right_mend_id:
        read_pair = reads_type1_right[read_pair_id]
        read = read_pair[0]
        read_mate = read_pair[1]
        new_seq = mend_read_part(ref, read, start, end, "type1", "right", svtype, subPos)
        qual = read.query_qualities
        read.query_sequence = new_seq
        read.query_qualities = qual
        total_reads.append([read, read_mate])
    print "deal type1 end......"
    return total_reads


def deal_type2(ref, reads_type2_left, reads_type2_right, freq, start, end, svtype, subPos=None):
    print "deal type2 start......"
    reads_left_num = len(reads_type2_left)
    reads_right_num = len(reads_type2_right)
    reads_left_mend_id = random_mendIDList(reads_left_num, freq)
    reads_right_mend_id = random_mendIDList(reads_right_num, freq)
    total_reads = []
    if svtype == "del":
        for read_pair_id in reads_left_mend_id:
            total_reads.append(reads_type2_left[read_pair_id])

        for read_pair_id in reads_right_mend_id:
            total_reads.append(reads_type2_right[read_pair_id])

    elif svtype == "inv":
        for read_pair_id in reads_left_mend_id:
            read_pair = reads_type2_left[read_pair_id]
            read = read_pair[0]
            read_mate = read_pair[1]
            new_seq = mend_read_part(ref, read, start, end, "type2", "left", svtype)
            qual = read.query_qualities
            read.query_sequence = new_seq
            read.query_qualities = qual
            total_reads.append([read, read_mate])

        for read_pair_id in reads_right_mend_id:
            read_pair = reads_type2_right[read_pair_id]
            read = read_pair[1]
            read_mate = read_pair[0]
            new_seq = mend_read_part(ref, read, start, end, "type2", "right", svtype)
            qual = read.query_qualities
            read.query_sequence = new_seq
            read.query_qualities = qual
            total_reads.append([read_mate, read])

    elif svtype == "dup":
        for read_pair_id in reads_left_mend_id:
            read_pair = reads_type2_left[read_pair_id]
            read = read_pair[0]
            new_seq = mend_read_part(ref, read, start, end, "type2", "left", svtype)
            qual = read.query_qualities
            read.query_sequence = new_seq
            read.query_qualities = qual
            read2 = read_pair[1]
            total_reads.append([read, read2])

        for read_pair_id in reads_right_mend_id:
            read_pair = reads_type2_right[read_pair_id]
            read2 = read_pair[1]
            new_seq = mend_read_part(ref, read2, start, end, "type2", "right", svtype)
            qual = read2.query_qualities
            read2.query_sequence = new_seq
            read2.query_qualities = qual
            read = read_pair[0]
            total_reads.append([read, read2])

    elif svtype == "trans_unbalance":
        for read_pair_id in reads_right_mend_id:
            read_pair = reads_type2_right[read_pair_id]
            read = read_pair[0]
            new_seq = mend_read_part(ref, read, start, end, "type2", "right", svtype, subPos=subPos)
            qual = read.query_qualities
            read.query_sequence = new_seq
            read.query_qualities = qual
            read2 = read_pair[1]
            total_reads.append([read, read2])

    elif svtype == "trans_balance":
        for read_pair_id in reads_left_mend_id:
            read_pair = reads_type2_left[read_pair_id]
            read = read_pair[0]
            new_seq = mend_read_part(ref, read, start, end, "type2", "left", svtype, subPos=subPos)
            qual = read.query_qualities
            read.query_sequence = new_seq
            read.query_qualities = qual
            read2 = read_pair[1]
            total_reads.append([read, read2])

        for read_pair_id in reads_right_mend_id:
            read_pair = reads_type2_right[read_pair_id]
            read2 = read_pair[1]
            new_seq = mend_read_part(ref, read2, start, end, "type2", "right", svtype, subPos=subPos)
            qual = read2.query_qualities
            read2.query_sequence = new_seq
            read2.query_qualities = qual
            read = read_pair[0]
            total_reads.append([read, read2])

    elif svtype == "trans_chrom":
        for read_pair_id in reads_left_mend_id:
            read_pair = reads_type2_left[read_pair_id]
            read = read_pair[0]
            new_seq = mend_read_part(ref, read, start, end, "type2", "left", svtype, subPos=subPos)
            qual = read.query_qualities
            read.query_sequence = new_seq
            read.query_qualities = qual
            read2 = read_pair[1]
            total_reads.append([read, read2])
    print "deal type2 end......"
    return total_reads


def deal_type3(reads_type3_left, reads_type3_right, freq, insertSize, start, end, svtype, supple1=None, supple2=None,
               subPos=None):
    # choose left read of start and right of end
    print "deal type3 start......"
    reads_left_num = len(reads_type3_left)
    reads_right_num = len(reads_type3_right)
    # print reads_left_num, reads_right_num
    reads_left_mend_id = random_mendIDList(reads_left_num, freq)
    reads_right_mend_id = random_mendIDList(reads_right_num, freq)
    total_del_reads = []
    total_add_reads = []
    total_modify_reads = []
    if svtype == "del":
        if len(supple1) == 0:
            print "Step3: Warning! No corresponding reads to pair"
            return [], [], []
        for read_pair_id in reads_left_mend_id:
            read_pair = reads_type3_left[read_pair_id]
            read = read_pair[0]
            try_time = 0
            while True:
                read_pair_tmp = random.sample(supple1, 1)[0]
                read_right = read_pair_tmp[1]
                insertSize_tmp = read_right.reference_end - read.reference_start - (end - start + 1)
                # print insertSize, insertSize_tmp
                if insertSize[0] <= insertSize_tmp <= insertSize[1] and check_reads_pair(read, read_right):
                    new_read = read_right
                    total_modify_reads.append([read, new_read])
                    total_del_reads.append([read_pair_tmp[0], read_pair[1]])
                    break
                try_time += 1
                if try_time > try_max_time:
                    print "can't find a mate read to match!"
                    break

    if svtype == "dup":
        if len(supple1) == 0:
            print "Step3: Warning! No corresponding reads to pair"
            return [], [], []
        for read_pair_id in reads_left_mend_id:
            read_pair = reads_type3_left[read_pair_id]
            read = read_pair[1]
            try_time = 0
            while True:
                read_left = random.sample(supple1, 1)[0][0]
                insertSize_tmp = end - read_left.reference_start + read.reference_end - start + 1
                if insertSize[0] <= insertSize_tmp <= insertSize[1] and check_reads_pair(read_left, read):
                    new_read = read_left
                    total_add_reads.append([read, new_read])
                    break
                try_time += 1
                if try_time > try_max_time:
                    print "can't find a mate read to match!"
                    break

    elif svtype == "inv":
        reads_type4 = supple1
        if len(reads_type4) == 0:
            print "Step3: Warning! No corresponding reads to pair"
            return [], [], []
        for read_pair_id in reads_left_mend_id:
            read_pair = reads_type3_left[read_pair_id]
            read = read_pair[0]
            try_time = 0
            while True:
                read_pair_tmp = copy.deepcopy(random.sample(reads_type4, 1)[0])
                read_right = read_pair_tmp[1]
                insertSize_tmp = start - read.reference_start + end - read_right.reference_start + 1
                if insertSize[0] <= insertSize_tmp <= insertSize[1] and check_reads_pair(read,
                                                                                         read_right):  ### reverse information?!!!
                    new_read = read_right
                    qual = new_read.query_qualities
                    new_read.query_sequence = getComplementarySeq(new_read.query_sequence)[::-1]
                    new_read.query_qualities = qual[::-1]
                    total_modify_reads.append([read, new_read])
                    total_modify_reads.append([read_pair_tmp[0], read_pair[1]])
                    break
                try_time += 1
                if try_time > try_max_time:
                    print "can't find a mate read to match!"
                    break

        for read_pair_id in reads_right_mend_id:
            read_pair = reads_type3_right[read_pair_id]
            read = read_pair[1]
            try_time = 0
            while True:
                read_pair_tmp = copy.deepcopy(random.sample(reads_type4, 1)[0])
                read_left = read_pair_tmp[0]
                insertSize_tmp = read_left.reference_end - start + read.reference_end - end + 1
                if insertSize[0] <= insertSize_tmp <= insertSize[1] and check_reads_pair(read_left, read):
                    new_read = read_left
                    qual = new_read.query_qualities
                    new_read.query_sequence = getComplementarySeq(new_read.query_sequence)[::-1]
                    new_read.query_qualities = qual[::-1]
                    total_modify_reads.append([new_read, read])
                    total_modify_reads.append([read_pair[0], read_pair_tmp[1]])
                    break
                try_time += 1
                if try_time > try_max_time:
                    print "can't find a mate read to match!"
                    break

    elif svtype == "trans_balance":
        reads_left_sub, reads_right_sub = supple1, supple2
        res = re.match("(\w*):(\d*)-(\d*)", subPos)
        trans_chr, trans_start, trans_end = res.group(1), int(res.group(2)), int(res.group(3))
        if len(reads_left_sub) == 0 or len(reads_right_sub) == 0:
            print "Step3: Warning! No corresponding reads to pair"
            return [], [], []
        for read_pair_id in reads_left_mend_id:
            read_pair = reads_type3_left[read_pair_id]
            read = read_pair[0]
            try_time = 0
            while True:
                read_pair_tmp = random.sample(reads_left_sub, 1)[0]
                read_right = read_pair_tmp[1]
                insertSize_tmp = start - read.reference_start + read_right.reference_end - trans_start + 1
                if insertSize[0] <= insertSize_tmp <= insertSize[1] and check_reads_pair(read, read_right):
                    new_read = read_right
                    total_modify_reads.append([read, new_read])
                    total_modify_reads.append([read_pair_tmp[0], read_pair[1]])
                    break
                try_time += 1
                if try_time > try_max_time:
                    print "can't find a mate read to match!"
                    break

        for read_pair_id in reads_right_mend_id:
            read_pair = reads_type3_right[read_pair_id]
            read = read_pair[1]
            try_time = 0
            while True:
                read_pair_tmp = random.sample(reads_right_sub, 1)[0]
                read_left = read_pair_tmp[0]
                insertSize_tmp = trans_end - read_left.reference_start + read.reference_end - end + 1
                if insertSize[0] <= insertSize_tmp <= insertSize[1] and check_reads_pair(read_left, read):
                    new_read = read_left
                    total_modify_reads.append([new_read, read])
                    total_modify_reads.append([read_pair[0], read_pair_tmp[1]])
                    break
                try_time += 1
                if try_time > try_max_time:
                    print "can't find a mate read to match!"
                    break
    elif svtype == "trans_chrom":
        reads_left_sub, reads_right_sub = supple1, supple2
        res = re.match("(\w*):(\d*)-(\d*)", subPos)
        trans_chr, trans_start, trans_end = res.group(1), int(res.group(2)), int(res.group(3))
        if len(reads_left_sub) == 0:
            print "Step3: Warning! No corresponding reads to pair"
            return [], [], []
        for read_pair_id in reads_left_mend_id:
            read_pair = reads_type3_left[read_pair_id]
            read = read_pair[0]
            try_time = 0
            while True:
                read_pair_tmp = random.sample(reads_left_sub, 1)[0]
                read_right = read_pair_tmp[1]
                insertSize_tmp = start - read.reference_start + read_right.reference_end - trans_start + 1
                if insertSize[0] <= insertSize_tmp <= insertSize[1] and check_reads_pair(read, read_right):
                    new_read = read_right
                    total_modify_reads.append([read, new_read])
                    total_modify_reads.append([read_pair_tmp[0], read_pair[1]])
                    break
                try_time += 1
                if try_time > try_max_time:
                    print "can't find a mate read to match!"
                    break
    elif svtype == "trans_unbalance":
        reads_left_sub, reads_right_sub = supple1, supple2
        res = re.match("(\w*):(\d*)-(\d*)", subPos)
        trans_chr, trans_start, trans_end = res.group(1), int(res.group(2)), int(res.group(3))
        if len(reads_left_sub) == 0 or len(reads_right_sub) == 0:
            print "Step3: Warning! No corresponding reads to pair"
            return [], [], []
        for read_pair_id in reads_left_mend_id:
            read_pair = reads_type3_left[read_pair_id]
            read = read_pair[0]
            try_time = 0
            while True:
                read_pair_tmp = random.sample(reads_left_sub, 1)[0]
                read_right = read_pair_tmp[1]
                insertSize_tmp = start - read.reference_start + read_right.reference_end - trans_start + 1
                if insertSize[0] <= insertSize_tmp <= insertSize[1] and check_reads_pair(read,
                                                                                         read_right):  # reverse information?!!!
                    new_read = read_right
                    total_modify_reads.append([read, new_read])
                    total_modify_reads.append([read_pair_tmp[0], read_pair[1]])
                    break
                try_time += 1
                if try_time > try_max_time:
                    print "can't find a mate read to match!"
                    break

        for read_pair_id in reads_right_mend_id:  # use left again, because right is none
            read_pair = reads_type3_right[read_pair_id]
            read = read_pair[1]
            try_time = 0
            while True:
                read_pair_tmp = random.sample(reads_right_sub, 1)[0]
                read_left = read_pair_tmp[0]
                insertSize_tmp = trans_end - read_left.reference_start + read.reference_end - end + 1
                if insertSize[0] <= insertSize_tmp <= insertSize[1] and check_reads_pair(read_left, read):
                    new_read = read_left
                    total_modify_reads.append([new_read, read])
                    total_modify_reads.append([read_pair[0], read_pair_tmp[1]])
                    break
                try_time += 1
                if try_time > try_max_time:
                    print "can't find a mate read to match!"
                    break
    print "deal type3 end......"
    return total_modify_reads, total_add_reads, total_del_reads


def deal_type4(reads_type4, freq, svtype, insertSize=None, cnvType=None):
    print "deal type4 start......"
    reads_num = len(reads_type4)
    reads_mend_id = random_mendIDList(reads_num, freq)
    total_reads = []
    if svtype == "del" or (svtype == "cnv" and cnvType == "loss"):
        for read_pair_id in reads_mend_id:
            total_reads.append(reads_type4[read_pair_id])

    elif svtype == "inv":
        total_reads = []

    elif svtype == "trans_balance":
        total_reads = []

    elif svtype == "dup" or (svtype == "cnv" and cnvType == "gain") or svtype == "trans_unbalance":
        reads_mend_left_id = random_mendIDList(reads_num / 2, freq)
        print "left mend count", len(reads_mend_left_id)
        for read_pair_id in reads_mend_left_id:
            # print read_pair_id
            read_pair = reads_type4[read_pair_id]
            read = read_pair[0]
            # print read.reference_start
            try_time = 0
            while True:
                read_right_pair_id = random.randint(read_pair_id + 1, min(read_pair_id + shift_num, reads_num - 1))
                read_right = reads_type4[read_right_pair_id][1]
                if read_right.query_name == read.query_name:
                    continue
                insertSize_tmp = read_right.reference_end - read.reference_start + 1
                if insertSize[0] <= insertSize_tmp <= insertSize[1] and check_reads_pair(read, read_right):
                    new_read = read_right
                    total_reads.append([read, new_read])
                    break
                try_time += 1
                if try_time > try_max_time:
                    print "can't find a mate read to match!"
                    break
                    # print "out loop1"

        reads_mend_right_id = [reads_num / 2 + i for i in random_mendIDList(reads_num - reads_num / 2, freq)]
        print "right mend count", len(reads_mend_right_id)
        for read_pair_id in reads_mend_right_id:
            # print read_pair_id
            read_pair = reads_type4[read_pair_id]
            read = read_pair[1]
            try_time = 0
            while True:
                read_left_pair_id = random.randint(max(read_pair_id - shift_num, 0), read_pair_id - 1)
                read_left = reads_type4[read_left_pair_id][0]
                if read_left.query_name == read.query_name:
                    continue
                insertSize_tmp = read.reference_end - read_left.reference_start + 1
                if insertSize[0] <= insertSize_tmp <= insertSize[1] and check_reads_pair(read_left, read):
                    new_read = read_left
                    total_reads.append([new_read, read])
                    break
                try_time += 1
                if try_time > try_max_time:
                    print "can't find a mate read to match!"
                    break
    print "deal type4 end......"
    return total_reads


def deal_type5(ref, reads_type5_left, reads_type5_right, freq, start, end, svtype, supple1=None, supple2=None,
               subPos=None):
    print "deal type5 start......"
    reads_left_num = len(reads_type5_left)
    reads_right_num = len(reads_type5_right)

    reads_left_mend_id = random_mendIDList(reads_left_num, freq)
    reads_right_mend_id = random_mendIDList(reads_right_num, freq)
    total_reads = []
    if svtype == "del":
        for read_pair_id in reads_left_mend_id:
            read_pair = reads_type5_left[read_pair_id]
            read_left = read_pair[0]
            read_right = read_pair[1]
            new_seq_left = mend_read_part(ref, read_left, start, end, "type5", "left", svtype)
            new_seq_right = mend_read_part(ref, read_right, start, end, "type5", "left", svtype)
            qual_left = read_left.query_qualities
            qual_right = read_right.query_qualities
            read_left.query_sequence = new_seq_left
            read_right.query_sequence = new_seq_right
            read_left.query_qualities = qual_left
            read_right.query_qualities = qual_right
            total_reads.append([read_left, read_right])

    elif svtype == "dup" or svtype == "trans_balance" or svtype == "inv" or svtype == "trans_unbalance":
        for read_pair_id in reads_left_mend_id:
            read_pair = reads_type5_left[read_pair_id]
            read_left = copy.deepcopy(read_pair[0])
            read_right = copy.deepcopy(read_pair[1])
            new_seq_left = mend_read_part(ref, read_left, start, end, "type5", "left", svtype, subPos=subPos)
            new_seq_right = mend_read_part(ref, read_right, start, end, "type5", "left", svtype, subPos=subPos)
            qual_left = read_left.query_qualities
            qual_right = read_right.query_qualities
            read_left.query_sequence = new_seq_left
            read_right.query_sequence = new_seq_right
            read_left.query_qualities = qual_left
            read_right.query_qualities = qual_right
            total_reads.append([read_left, read_right])
        for read_pair_id in reads_right_mend_id:
            read_pair = reads_type5_right[read_pair_id]
            read_left = copy.deepcopy(read_pair[0])
            read_right = copy.deepcopy(read_pair[1])
            new_seq_left = mend_read_part(ref, read_left, start, end, "type5", "right", svtype, subPos=subPos)
            new_seq_right = mend_read_part(ref, read_right, start, end, "type5", "right", svtype, subPos=subPos)
            qual_left = read_left.query_qualities
            qual_right = read_right.query_qualities
            read_left.query_sequence = new_seq_left
            read_right.query_sequence = new_seq_right
            read_left.query_qualities = qual_left
            read_right.query_qualities = qual_right
            total_reads.append([read_left, read_right])
    elif svtype == "trans_chrom":
        for read_pair_id in reads_left_mend_id:
            read_pair = reads_type5_left[read_pair_id]
            read_left = copy.deepcopy(read_pair[0])
            read_right = copy.deepcopy(read_pair[1])
            new_seq_left = mend_read_part(ref, read_left, start, end, "type5", "left", svtype, subPos=subPos)
            new_seq_right = mend_read_part(ref, read_right, start, end, "type5", "left", svtype, subPos=subPos)
            qual_left = read_left.query_qualities
            qual_right = read_right.query_qualities
            read_left.query_sequence = new_seq_left
            read_right.query_sequence = new_seq_right
            read_left.query_qualities = qual_left
            read_right.query_qualities = qual_right
            total_reads.append([read_left, read_right])

    print "deal type5 end......"

    return total_reads


def random_mendIDList(totalReadsNum, frac):
    if totalReadsNum == 0:
        return []
    if frac <= 1:
        mendList = random.sample(range(totalReadsNum), int(totalReadsNum * frac))
        return mendList
    elif frac > 1:
        cnt = int(frac)
        frac_sub = frac / (cnt + 1)
        total_mendList = []
        for i in range(cnt + 1):
            mendList = random.sample(range(totalReadsNum), int(totalReadsNum * frac_sub))
            total_mendList.extend(mendList)
        return total_mendList


def deal_type6(ref, reads_type1_left, reads_type1_right, freq, start, end, svtype, subPos=None):
    print "deal type6 start......"
    reads_left_num = len(reads_type1_left)
    reads_right_num = len(reads_type1_right)
    reads_left_mend_id = random_mendIDList(reads_left_num, freq)
    reads_right_mend_id = random_mendIDList(reads_right_num, freq)

    total_reads = []
    for read_pair_id in reads_left_mend_id:
        read_pair = reads_type1_left[read_pair_id]
        read = read_pair[0]
        new_seq = None
        if svtype == "del" or svtype == "inv" or svtype == "trans_balance" or svtype == "trans_unbalance":
            if start - read.reference_start >= read.reference_end - start:
                new_seq = mend_read_part(ref, read, start, end, "type6", "left", svtype, subPos=subPos)
            else:
                continue
        elif svtype == "dup":
            if start - read.reference_start <= read.reference_end - start:
                new_seq = mend_read_part(ref, read, start, end, "type6", "left", svtype, subPos=subPos)
            else:
                continue
        if new_seq:
            qual = read.query_qualities
            read.query_sequence = new_seq
            read.query_qualities = qual
        total_reads.append([read])

    for read_pair_id in reads_right_mend_id:
        read_pair = reads_type1_right[read_pair_id]
        read = read_pair[0]
        new_seq = None
        if svtype == "del" or svtype == "inv" or svtype == "trans_balance" or svtype == "trans_unbalance":
            if end - read.reference_start <= read.reference_end - end:
                new_seq = mend_read_part(ref, read, start, end, "type6", "right", svtype, subPos=subPos)
            else:
                continue
        elif svtype == "dup":
            if end - read.reference_start >= read.reference_end - end:
                new_seq = mend_read_part(ref, read, start, end, "type6", "right", svtype, subPos=subPos)
            else:
                continue
        if new_seq:
            qual = read.query_qualities
            read.query_sequence = new_seq
            read.query_qualities = qual
        total_reads.append([read])
    print "deal type6 end......"
    return total_reads


def deal_type7(reads_type4, freq, svtype):
    print "deal type7 start......"
    reads_num = len(reads_type4)
    reads_mend_id = random_mendIDList(reads_num, freq)
    total_reads = []
    if svtype == "del" or "dup" or "trans_unbalance" or svtype == "cnv":
        for read_pair_id in reads_mend_id:
            total_reads.append(reads_type4[read_pair_id])

    return total_reads
