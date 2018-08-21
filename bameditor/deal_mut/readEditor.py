import random
import pysam
from array import array


def editRead(read, refFile, mutList):
    seq = read.query_sequence
    qual = read.query_qualities
    if qual is None:
        print "Warning: " + str(read) + " quality is None. Continue!"
        return False
    ref = pysam.FastaFile(refFile)
    seq_len = len(seq)
    shift = 0
    for mut in mutList:
        res = None
        if mut.muttype == "snv":
            res = snvEdit(read, mut.start)
            if res:
                pos = res
                seq = seq[:pos + shift] + mut.alt + seq[pos + 1 + shift:]
        elif mut.muttype == "ins":
            res = insEdit(read, mut.start)
            if res:
                pos = res
                seq = seq[:pos + 1 + shift] + mut.alt + seq[pos + 1 + shift:]
                # print qual, pos+1+shift-2,pos+1+shift+2
                pos_s = pos + 1 + shift - 2
                pos_e = pos + 1 + shift + 2
                if pos_s < 0:
                    pos_s = 0
                if len(qual[pos_s:pos_e]) != 0:
                    max_qual, min_qual = max(qual[pos_s:pos_e]), min(qual[pos_s:pos_e])
                else:
                    max_qual, min_qual = max(qual), min(qual)
                ins_qual = []
                for i in range(len(mut.alt)):
                    ins_qual.append(random.randint(min_qual, max_qual))
                ins_qual = array('B', ins_qual)
                qual = qual[:pos + 1 + shift] + ins_qual + qual[pos + 1 + shift:]
                shift = shift + len(mut.alt)
        elif mut.muttype == "del":
            res = delEdit(read, mut.start, mut.end)
            if res:
                start, end = res
                seq = seq[:start + shift] + seq[end + 1 + shift:]
                qual = qual[:start + shift] + qual[end + 1 + shift:]
                shift = shift - (end - start + 1)
        elif mut.muttype == "sub":
            res = subEdit(read, mut.start, mut.end)
            if res:
                start, end = res
                # print "start.........................."
                seq = seq[:start + shift] + mut.alt + seq[end + 1 + shift:]
                # print seq
                # print "end............................"
                if len(qual[start + shift:end + 1 + shift]) != 0:
                    max_qual, min_qual = max(qual[start + shift:end + 1 + shift]), min(
                        qual[start + shift:end + 1 + shift])
                else:
                    max_qual, min_qual = max(qual), min(qual)
                ins_qual = []
                for i in range(len(mut.alt)):
                    ins_qual.append(random.randint(min_qual, max_qual))
                ins_qual = array('B', ins_qual)
                qual = qual[:start + shift] + ins_qual + qual[end + 1 + shift:]
                shift = shift - (end - start + 1) + len(mut.alt)

        if not res:
            # print read, mut.muttype
            return False
    if shift >= 0:
        sequence = seq[:seq_len]
        quality = qual[:seq_len]
    else:
        start = read.reference_end
        end = start - shift
        fill_seq = ref.fetch(read.reference_name, start, end)
        sequence = seq + fill_seq
        if len(qual[-6:]) != 0:
            max_qual, min_qual = max(qual[-6:]), min(qual[-6:])
        else:
            max_qual, min_qual = max(qual), min(qual)
        fill_qual = array('B')
        for i in range(-shift):
            fill_qual.append(random.randint(min_qual, max_qual))
        quality = qual + fill_qual
    return sequence, quality, shift


def snvEdit(read, pos):
    posPairList = read.get_aligned_pairs
    if pos < read.reference_start or pos > read.reference_end:
        print "Warning: SNV is out of this read!"
        return False

    mutPosIndex = None
    for i, pair in enumerate(posPairList):
        if pair[1] == pos:
            mutPosIndex = pair[0]
            break
    if mutPosIndex is None:
        print "Warning: A deletion existed of this read: %s at pos %s, can't spike in snv" % (read.query_name, pos)
        return False
    return mutPosIndex


def insEdit(read, pos):
    posPairList = read.get_aligned_pairs
    if pos + 1 < read.reference_start or pos > read.reference_end:
        print "Warning: ins is out of this read!"
        return False
    elif pos + 1 == read.reference_start:
        mutPosIndex = -1
    else:
        mutPosIndex = None
        for i, pair in enumerate(posPairList):
            if pair[1] == pos:
                mutPosIndex = pair[0]
        if mutPosIndex is None:
            print "Warning: A deletion existed of this read: %s at pos %s, can't spike in ins" % (read.query_name, pos)
            return False
    return mutPosIndex


def delEdit(read, start, end):
    posPairList = read.get_aligned_pairs
    if (start < read.reference_start and end < read.reference_start) or (
            start > read.reference_end and end > read.reference_end):
        print "Warning: del is out of this read!"
        return False

    elif start < read.reference_start <= end:
        start = read.reference_start

    elif start <= read.reference_end < end:
        end = read.reference_end

    startPosIndex, endPosIndex = None, None
    for pair in posPairList:
        if pair[1] <= start:
            if pair[0] is not None:
                startPosIndex = pair[0]
        if pair[1] >= start and pair[1] <= end:
            if pair[0] is not None:
                endPosIndex = pair[0]
        if pair[1] > end:
            break
    if startPosIndex is None or endPosIndex is None:
        print "Warning: A deletion existed of this read: %s at pos %s-%s, can't spike in del" % (
        read.query_name, startPosIndex, endPosIndex)
        return False

    return (startPosIndex, endPosIndex)


def subEdit(read, start, end):
    posPairList = read.get_aligned_pairs
    if (start < read.reference_start and end < read.reference_start) or (
            start > read.reference_end and end > read.reference_end):
        print "Warning: sub is out of this read!"
        return False

    elif start < read.reference_start <= end:
        start = read.reference_start

    elif start <= read.reference_end < end:
        end = read.reference_end

    startPosIndex, endPosIndex = None, None
    # print start, end
    for pair in posPairList:
        if pair[1] <= start:
            if pair[0] is not None:
                startPosIndex = pair[0]
        if pair[1] >= start and pair[1] <= end:
            if pair[0] is not None:
                endPosIndex = pair[0]
        if pair[1] > end:
            break
    if startPosIndex is None or endPosIndex is None:
        print "Warning: A deletion existed of this read: %s at pos %s-%s, can't spike in sub" % (
        read.query_name, startPosIndex, endPosIndex)
        return False
    # print startPosIndex, endPosIndex
    return (startPosIndex, endPosIndex)
