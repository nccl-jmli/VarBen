# Author: Shuangsang Fang
# Date: 2018/8/21
import sys
from uuid import uuid4

import math
import random
from array import array
import pysam


def count_coverage(bam, chr, start, end):
    coverage = 0
    for read in bam.fetch(chr, start, end + 1):
        coverage += 1
    return coverage


def calPolymorphism(bam, chrom, start, end):
    reads = bam.fetch(chrom, start, end + 1)
    numDict = {}
    readsNum = 0
    for pos in range(start, end):
        numDict[pos] = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    for read in reads:
        readsNum += 1
        for pair in read.get_aligned_pairs():
            if start <= pair[1] < end:
                pos = pair[0]
                if pos is None:
                    continue
                if pos < read.query_length:
                    nul = read.query_sequence[pos]
                    if nul in 'ATCG':
                        numDict[pair[1]][nul] += 1
    return numDict, readsNum


def getReadStrand(read):
    if read.is_read1:
        return '1'
    elif read.is_read2:
        return '2'
    else:
        return '0'


def getKeyName(read):
    return read.query_name + "/" + getReadStrand(read)


def getMateKeyName(keyName):
    tmp = keyName.split("/")
    if tmp[-1] == '1':
        return "".join(tmp[:-1]) + "/2"
    elif tmp[-1] == '2':
        return "".join(tmp[:-1]) + "/1"
    else:
        return "".join(tmp[:-1]) + "/0"


def find_mate(read, bam):
    """ AlignmentFile.mate() can return a non-primary alignment, so use this function instead """
    chrom = read.next_reference_name
    for rec in bam.fetch(chrom, read.next_reference_start, read.next_reference_start + 1):
        if rec.query_name == read.query_name and rec.reference_start == read.next_reference_start:
            if not rec.is_secondary and bin(rec.flag & 2048) != bin(2048):
                if rec.is_read1 != read.is_read1:
                    return rec
    return None


def check_read(read, is_single):
    if is_single:
        pass
    else:
        if read.flag in [163, 83, 99, 147]:
            return True
        else:
            return False
    return False


def check_reads_pair(read_left, read_right):
    if (not read_left.is_reverse and read_right.is_reverse and (
                (read_left.is_read2 and read_right.is_read1) or (read_left.is_read1 and read_right.is_read2))):
        return True
    else:
        return False


def get_new_readname():
    newId = str(uuid4())
    return newId


def get_insertSize_range(bam_file, readLength, is_single):
    return [100, 1000]
    # if is_single:
    #     return []
    # import numpy
    # total_insertSize = []
    # bam = pysam.AlignmentFile(bam_file, 'rb')
    # fout = open("insertSize.txt", 'w')
    # for read in bam.fetch():
    #     if read.cigarstring == str(readLength) + "M" and read.flag in [163, 83, 99, 147]:
    #         insertSize = abs(read.template_length)
    #         total_insertSize.append(insertSize)
    #         fout.write(str(insertSize) + "\n")
    # fout.close()
    # total_insertSize.sort()
    # total_len = len(total_insertSize)
    # total_insertSize_narray = numpy.array(total_insertSize)
    # sum1 = total_insertSize_narray.sum()
    # total_insertSize_narray2 = total_insertSize_narray * total_insertSize_narray
    # sum2 = total_insertSize_narray2.sum()
    # mean = sum1 / total_len
    # stdev = math.sqrt(sum2 / total_len - mean ** 2)
    #
    # min_insert = mean - stdev * 3 + 100 if mean - stdev * 3 > 0 else 100
    # max_insert = mean + stdev * 3 + 100 if mean - stdev * 3 > 400 else 500
    # # print mean, var
    # return [min_insert, max_insert]


def get_mate(read, bam):
    if read.is_paired and not read.mate_is_unmapped:
        return bam.mate(read)


def getComplementarySeq(seq):
    complementDict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    new_seq = ""
    for nul in seq:
        nul = nul.upper()
        if nul in complementDict:
            new_seq = new_seq + complementDict[nul]
        else:
            new_seq = new_seq + nul
    return new_seq


def reCalculateFlowSignal(sequence, origFlowSignal, flowOrder, libKey, barcode, reverse):
    seqPos = 0
    # atcg = 'ATCG'
    # randomSeq = ""
    # for i in range(7):
    #     randomID = random.randint(0, 3)
    #     randomSeq += atcg[randomID]

    randomSeq = "ATCACCG"
    sequence = sequence.upper()
    if reverse:
        rev_seq = getComplementarySeq(sequence[::-1])
        new_seq = libKey + barcode + rev_seq + randomSeq
    else:
        new_seq = libKey + barcode + sequence + randomSeq

    flowPos = 0
    max_signal = 0
    min_signal = 1000
    flowSize = 500
    newFlowSignal = array('h')
    while seqPos < len(new_seq):
        # flowValue = flowOrder[flowPos]
        if flowOrder[flowPos] == new_seq[seqPos]:
            if seqPos + 1 < len(libKey):
                if origFlowSignal[flowPos] > max_signal:
                    max_signal = origFlowSignal[flowPos]
                if origFlowSignal[flowPos] < min_signal:
                    min_signal = origFlowSignal[flowPos]
                newFlowSignal.append(origFlowSignal[flowPos])
                seqPos += 1
            # elif seqPos + 1 < len(libKey) + len(barcode):
            #     newFlowSignal.append(origFlowSignal[flowPos])
            #     seqPos += 1
            else:
                tmp_singal = 0
                while seqPos < len(new_seq) and new_seq[seqPos] == flowOrder[flowPos]:
                    tmp_singal += random.randint(min_signal, max_signal)
                    seqPos += 1
                newFlowSignal.append(tmp_singal)
        else:
            newFlowSignal.append(random.randint(-30, 30))
        if len(newFlowSignal) == flowSize:
            break
        flowPos += 1
        if flowPos == len(flowOrder):
            flowPos = 0

    start = len(libKey + barcode)
    if seqPos > len(libKey + barcode + sequence.upper()):
        seqPos = len(libKey + barcode + sequence.upper())

    final_seq = new_seq[start:seqPos]
    if reverse:
        final_seq = getComplementarySeq(final_seq[::-1])
    if len(newFlowSignal) < flowSize:
        newFlowSignal.append(random.randint(-30, 30))
    # else:
    #     final_seq = final_seq[:-7]
    # newFlowSignalStr = ",".join([singal[0]] + [str(x) for x in newFlowSignal])
    return newFlowSignal, final_seq


def reCalculateSequence(origFlowSignal, flowOrder, libKey, barcode):
    seq = ""
    i = 0
    pos = 0
    flowOrderN = flowOrder * 50
    while pos < len(origFlowSignal):
        if int(origFlowSignal[pos]) >= 190:
            n = int(origFlowSignal[pos]) / 190
            seq += flowOrderN[i] * n
            # print flowOrderN[i]

        i += 1
        pos += 1

    print seq, len(libKey + barcode), len(libKey), len(barcode)
    print seq[len(libKey + barcode):]


def deal_life_reads(read, flow_order, lib_key, barcode):
    """deal life reads, modify the flow signal
    :param read: read object
    :param flow_order: flow order
    :param lib_key: lib key
    :param barcode: barcode
    :return: modified read
    """
    orig_flow_signal = None
    for tag in read.tags:
        if tag[0] == "ZM":
            orig_flow_signal = tag[1]
            break
    if not orig_flow_signal:
        print "Please make sure that it is a life read, since there is no tag ZM"
        return False

    flow_signal_str, final_seq = reCalculateFlowSignal(read.query_sequence, orig_flow_signal, flow_order,
                                                       lib_key, barcode, read.is_reverse)

    # change tags! change ZM
    tags = read.tags
    new_tags = []
    for tag in tags:
        if tag[0] == "ZM":
            new_tags.append((tag[0], flow_signal_str))
        else:
            new_tags.append(tag)
    read.tags = new_tags
    qual = read.query_qualities
    read.query_sequence = final_seq
    read.query_qualities = qual[:len(final_seq)]
    return read


def add_tag(read):
    tags = read.tags
    tags.append(("BE", 1))
    read.tags = tags
    return read
