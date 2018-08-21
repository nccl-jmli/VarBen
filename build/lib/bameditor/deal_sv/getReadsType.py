from multiprocessing.pool import Pool
import pysam
from bameditor.common import bamobject
from bameditor.common.bamobject import Read
from bameditor.common.methods import find_mate


def read_bam_reads(bamFileList):
    total_reads_pair = []
    for bamfile in bamFileList:
        bam = pysam.AlignmentFile(bamfile, 'rb')
        read_left, read_right = None, None
        for read in bam.fetch():
            if read_left == None:
                read_left = read
            elif read_right == None:
                read_right = read
                total_reads_pair.append([read_left, read_right])
                read_left, read_right = read, None
    return total_reads_pair


def getPosType_pair_pos(read1, read2, pos):
    # 0 represent no intersection, 1 represent left read intersection, 2 represent right intersection, 3 represent no intersection point but overlap
    # -1 represent reads all in the left of pos
    # -2 represent reads all in the right of pos
    # if read1.reference_start > read2.reference_start:
    #     read1, read2 = read2, read1
    read1_start, read1_end = read1.reference_start, read1.reference_end
    read2_start, read2_end = read2.reference_start, read2.reference_end
    if read1_end < pos:
        if read2_end < pos:
            return -1
        elif read2_start <= pos <= read2_end:
            return 2
        elif read2_start > pos:
            return 3
    elif read1_start <= pos <= read1_end:
        if read2_start <= pos <= read2_end:
            return 4
        else:
            return 1
    else:
        if read1_start > pos:
            return -2


def getPosType_pair(read1, read2, start, end):
    type_start = getPosType_pair_pos(read1, read2, start)
    if type_start > 0:
        return "left", type_start

    else:
        type_end = getPosType_pair_pos(read1, read2, end)

        if type_end > 0:
            return "right", type_end
        else:
            if type_start == -2 and type_end == -1:
                return "center", 0
            elif (type_start == -1 and type_end == -1) or (type_start == -2 and type_end == -2):
                return "side", 0


def getPosType_single_pos(read1, pos):
    # 0 represent read intersection
    # -1 represent reads all in the left of pos
    # -2 represent reads all in the right of pos
    read1_start = read1.reference_start
    read1_end = read1.reference_end
    if read1_end < pos:
        return -1
    elif read1_start <= pos <= read1_end:
        return 1
    else:
        return -2


def getPosType_single(read1, start, end):
    type_start = getPosType_single_pos(read1, start)
    if type_start > 0:
        return "left", type_start
    else:
        type_end = getPosType_single_pos(read1, end)
        if type_end > 0:
            return "right", type_end
        else:
            if type_start == -1 or type_end == -2:
                return "side", 0
            elif type_start == -2 and type_end == -1:
                return "center", 0
            else:
                raise Exception("Impossible!")


def posType_sub_single(bamFile, chr, sub_start, sub_end, start, end, minmapq, is_multmapfilter):
    bam = pysam.AlignmentFile(bamFile, 'rb')
    reads = bam.fetch(chr, sub_start, sub_end)

    total_reads = []
    for read in reads:
        total_reads.append(read)
    print len(total_reads)
    reads_type6_left = []  # 6. in left place of del and read is on the breakpoint
    reads_type6_right = []  # 6. in right place of del and read is on the breakpoint
    reads_type7 = []  # 7. reads within the del
    filtered_reads_num = 0

    for read in total_reads:
        if read.is_secondary or int(read.mapping_quality) < minmapq:
            filtered_reads_num += 1
            continue
        position, pos_type = getPosType_single(read, start, end)
        read = Read(read)
        if position == "left":
            if pos_type == 1:
                reads_type6_left.append([read])
            else:
                raise Exception("Impossible!")
        elif position == "right":
            if pos_type == 1:
                reads_type6_right.append([read])
            else:
                raise Exception("Impossible!")
        elif position == "center":
            reads_type7.append([read])

    bam.close()
    return reads_type6_left, reads_type6_right, reads_type7, filtered_reads_num


def posType_sub_paired(bamFile, chr, sub_start, sub_end, start, end, readLength, minmapq, is_multmapfilter,
                       extension=None):
    bam = pysam.AlignmentFile(bamFile, 'rb')
    print sub_start, sub_end
    reads = bam.fetch(chr, sub_start, sub_end)

    total_reads = []
    for read in reads:
        total_reads.append(read)
    reads_type1_left = []  # 1. in left place of del and second read is on the breakpoint
    reads_type1_right = []  # 1. in right place of del and first read is on the breakpoint
    reads_type2_left = []  # 2. in left place of del and first read is on the breakpoint
    reads_type2_right = []  # 2. in right place of del and second read is on the breakpoint
    reads_type3_left = []  # 3. left read and right read crossover start breakpoint with no intersection
    reads_type3_right = []  # 3. left read and right read crossover end breakpoint with no intersection
    reads_type4 = []  # 4. reads within the del
    reads_type5_left = []  # 5. in left place of del and first read and right read are all has intersection
    reads_type5_right = []  # 3. in right place of del and first read and right read are all has intersection
    filtered_reads_num = 0
    print "Sub reads num: ", len(total_reads)
    reverse_num = 0
    for read in total_reads:
        mate_read = find_mate(read, bam)
        if not (mate_read and read.cigarstring == str(readLength) + "M" and mate_read.cigarstring == str(
                readLength) + "M"
                and read.flag in [163, 83, 99, 147] and mate_read.flag in [163, 83, 99, 147]) \
                or (int(read.mapping_quality) < minmapq or int(
                    mate_read.mapping_quality) < minmapq):  # check_reads_pair(read, mate_read))
            filtered_reads_num += 1
            continue

        if is_multmapfilter and (read.has_tag("XA") or mate_read.has_tag("XA")):
            filtered_reads_num += 1
            continue

        if not extension:
            if read.is_reverse:
                if mate_read.reference_end >= start:
                    reverse_num += 1
                    continue
                else:
                    read, mate_read = mate_read, read
        else:
            if read.is_reverse:
                continue

        position, pos_type = getPosType_pair(read, mate_read, start, end)
        read_left = Read(read)
        read_right = Read(mate_read)
        if position == "left":
            if pos_type == 1:
                reads_type2_left.append([read_left, read_right])
            elif pos_type == 2:
                reads_type1_left.append([read_left, read_right])
            elif pos_type == 3:
                reads_type3_left.append([read_left, read_right])
            elif pos_type == 4:
                reads_type5_left.append([read_left, read_right])

        elif position == "right":
            if pos_type == 1:
                reads_type1_right.append([read_left, read_right])
            elif pos_type == 2:
                reads_type2_right.append([read_left, read_right])
            elif pos_type == 3:
                reads_type3_right.append([read_left, read_right])
            elif pos_type == 4:
                reads_type5_right.append([read_left, read_right])
        elif position == "center":
            reads_type4.append([read_left, read_right])
    print "Reverse num: ", reverse_num
    bam.close()
    return reads_type1_left, reads_type1_right, reads_type2_left, reads_type2_right, reads_type3_left, reads_type3_right, reads_type4, reads_type5_left, reads_type5_right, filtered_reads_num


def pos_type_classify(bamfile, chrom, start, end, is_single, read_length, temp_dir, extension=None, center=True,
                      maxsize=None, process=20, minmapq=0, is_multmapfilter=False):
    print bamfile, chrom, start, end, is_single, read_length, temp_dir, extension, center
    if is_single:
        total_reads_type6_left = []  # 6. in left place of del and second read is on the breakpoint
        total_reads_type6_right = []  # 6. in right place of del and first read is on the breakpoint
        total_reads_type7 = []  # 7. reads within the del
        # temp_prefix = "%s/classify_%s" % (temp_dir, sub_num)
        if extension:
            rel_start = start - extension
            rel_end = end + extension
        else:
            rel_start = start
            rel_end = end
        if center:
            reads_type6_left, reads_type6_right, reads_type7, filtered_reads_num = posType_sub_single(bamfile, chrom,
                                                                                                      rel_start,
                                                                                                      rel_end, start,
                                                                                                      end, minmapq,
                                                                                                      is_multmapfilter)
        else:
            rel_start_left = rel_start
            rel_end_left = start + maxsize
            rel_start_right = end - maxsize
            rel_end_right = rel_end
            reads_type6_left_1, reads_type6_right_1, reads_type7_1, filtered_reads_num_1 = posType_sub_single(bamfile,
                                                                                                              chrom,
                                                                                                              rel_start_left,
                                                                                                              rel_end_left,
                                                                                                              start,
                                                                                                              end,
                                                                                                              minmapq,
                                                                                                              is_multmapfilter)
            reads_type6_left_2, reads_type6_right_2, reads_type7_2, filtered_reads_num_2 = posType_sub_single(bamfile,
                                                                                                              chrom,
                                                                                                              rel_start_right,
                                                                                                              rel_end_right,
                                                                                                              start,
                                                                                                              end,
                                                                                                              minmapq,
                                                                                                              is_multmapfilter)
            reads_type6_left = reads_type6_left_1 + reads_type6_right_1
            reads_type6_right = reads_type6_right_1 + reads_type6_right_2
            reads_type7 = reads_type7_1 + reads_type7_2
            filtered_reads_num = filtered_reads_num_1 + filtered_reads_num_2

        total_reads_type6_left.extend(reads_type6_left)
        total_reads_type6_right.extend(reads_type6_right)
        total_reads_type7.extend(reads_type7)
        total_filtered_reads = filtered_reads_num
        print total_reads_type6_left, total_reads_type6_right, total_reads_type7, total_filtered_reads
        return total_reads_type6_left, total_reads_type6_right, total_reads_type7, total_filtered_reads
    else:
        total_reads_type1_left = []  # 1. in left place of del and second read is on the breakpoint
        total_reads_type1_right = []  # 1. in right place of del and first read is on the breakpoint
        total_reads_type2_left = []  # 2. in left place of del and first read is on the breakpoint
        total_reads_type2_right = []  # 2. in right place of del and second read is on the breakpoint
        total_reads_type3_left = []  # 3. in left place of del and first read and right read is crossover breakpoint with no intersection
        total_reads_type3_right = []  # 3. in right place of del and first read and right read is crossover breakpoint with no intersection
        total_reads_type4 = []  # 4. reads within the del
        total_reads_type5_left = []  # 5. in left place of del and first read and right read are all has intersection
        total_reads_type5_right = []  # 3. in right place of del and first read and right read are all has intersection
        total_filtered_reads = 0

        length = end - start + 1
        sub_num = length / read_length

        # when start = end, translocation of chromosome
        if start == end:
            rel_start = start - maxsize
            rel_end = end + maxsize
            print rel_start, rel_end
            # temp_prefix = "%s/classify_%s" % (temp_dir, "whole")
            (reads_type1_left, reads_type1_right, reads_type2_left, reads_type2_right, reads_type3_left,
             reads_type3_right, reads_type4, reads_type5_left, reads_type5_right, filtered_reads_num) \
                = posType_sub_paired(bamfile, chrom, rel_start, rel_end, start, end, read_length, minmapq,
                                     is_multmapfilter, extension=extension)
            total_reads_type1_left.extend(reads_type1_left)
            total_reads_type1_right.extend(reads_type1_right)
            total_reads_type2_left.extend(reads_type2_left)
            total_reads_type2_right.extend(reads_type2_right)
            total_reads_type3_left.extend(reads_type3_left)
            total_reads_type3_right.extend(reads_type3_right)
            total_reads_type4.extend(reads_type4)
            total_reads_type5_left.extend(reads_type5_left)
            total_reads_type5_right.extend(reads_type5_right)
            total_filtered_reads = filtered_reads_num

        # end - start < read_length and there is no need to extend its scope
        elif sub_num == 0 and not extension:
            # temp_prefix = "%s/classify_%s" % (temp_dir, sub_num)
            (reads_type1_left, reads_type1_right, reads_type2_left, reads_type2_right, reads_type3_left,
             reads_type3_right, reads_type4, reads_type5_left, reads_type5_right, filtered_reads_num) \
                = posType_sub_paired(bamfile, chrom, start, end, start, end, read_length, minmapq, is_multmapfilter,
                                     extension=extension)
            total_reads_type1_left.extend(reads_type1_left)
            total_reads_type1_right.extend(reads_type1_right)
            total_reads_type2_left.extend(reads_type2_left)
            total_reads_type2_right.extend(reads_type2_right)
            total_reads_type3_left.extend(reads_type3_left)
            total_reads_type3_right.extend(reads_type3_right)
            total_reads_type4.extend(reads_type4)
            total_reads_type5_left.extend(reads_type5_left)
            total_reads_type5_right.extend(reads_type5_right)
            total_filtered_reads = filtered_reads_num

        # there should be more than one process to calculate.
        else:
            run_pool = Pool(process)
            result_list = []
            # extension the range to cover whole reads
            if extension:
                rel_start = start - extension
                rel_end = end + extension
                length = rel_end - rel_start + 1
                sub_num = length / read_length
            else:
                rel_start = start
                rel_end = end
            # if center should be consider or center is no need to consider, but the center size is too less
            if center or (not center and maxsize is not None and length < maxsize * 2):
                for i in range(sub_num):
                    sub_start = i * read_length + rel_start
                    if i == sub_num - 1:
                        sub_end = rel_end
                    else:
                        sub_end = sub_start + 1
                    print "Sub Process: %s" % i, sub_start, sub_end
                    result_list.append(
                        run_pool.apply_async(posType_sub_paired, args=(
                            bamfile, chrom, sub_start, sub_end, start, end, read_length, minmapq, is_multmapfilter,
                            extension))
                    )
                run_pool.close()
                run_pool.join()
            # if center is no need to consider
            else:
                rel_start_left = rel_start
                rel_end_left = start + maxsize
                rel_start_right = end - maxsize
                rel_end_right = rel_end
                # print rel_start_left, rel_end_left, rel_start_right, rel_end_right
                length = rel_end_left - rel_start_left + 1
                sub_num = length / read_length
                for i in range(sub_num):
                    sub_start = i * read_length + rel_start_left
                    if i == sub_num - 1:
                        sub_end = rel_end_left
                    else:
                        sub_end = sub_start + 1
                    print "Sub Process: %s" % i, sub_start, sub_end
                    # temp_prefix = "%s/classify_%s" % (temp_dir, i)
                    result_list.append(
                        run_pool.apply_async(posType_sub_paired, args=(
                            bamfile, chrom, sub_start, sub_end, start, end, read_length, minmapq, is_multmapfilter,
                            extension))
                    )
                length = rel_end_right - rel_start_right + 1
                sub_num = length / read_length
                for i in range(sub_num):
                    sub_start = i * read_length + rel_start_right
                    if i == sub_num - 1:
                        sub_end = rel_end_right
                    else:
                        sub_end = sub_start + 1
                    print "Sub Process: %s" % i, sub_start, sub_end
                    # temp_prefix = "%s/classify_%s" % (temp_dir, i)
                    result_list.append(
                        run_pool.apply_async(posType_sub_paired, args=(
                            bamfile, chrom, sub_start, sub_end, start, end, read_length, minmapq, is_multmapfilter,
                            extension))
                    )
                run_pool.close()
                run_pool.join()

            for res in result_list:
                reads_type1_left, reads_type1_right, reads_type2_left, reads_type2_right, reads_type3_left, reads_type3_right, reads_type4, reads_type5_left, reads_type5_right, filtered_reads_num = res.get()
                total_reads_type1_left.extend(reads_type1_left)
                total_reads_type1_right.extend(reads_type1_right)
                total_reads_type2_left.extend(reads_type2_left)
                total_reads_type2_right.extend(reads_type2_right)
                total_reads_type3_left.extend(reads_type3_left)
                total_reads_type3_right.extend(reads_type3_right)
                total_reads_type4.extend(reads_type4)
                total_reads_type5_left.extend(reads_type5_left)
                total_reads_type5_right.extend(reads_type5_right)
                total_filtered_reads += filtered_reads_num

        print "type1_left: %s; type1_right: %s, type2_left: %s; type2_right: %s, type3_left: %s; " \
              "type3_right: %s, type4: %s; type5_left: %s; type5_right: %s" % (
                  len(total_reads_type1_left), len(total_reads_type1_right), len(total_reads_type2_left),
                  len(total_reads_type2_right), len(total_reads_type3_left), len(total_reads_type3_right),
                  len(total_reads_type4), len(total_reads_type5_left), len(total_reads_type5_right))
        print "total_filtered_reads: %s" % total_filtered_reads
        return total_reads_type1_left, total_reads_type1_right, total_reads_type2_left, total_reads_type2_right, total_reads_type3_left, total_reads_type3_right, total_reads_type4, total_reads_type5_left, total_reads_type5_right, total_filtered_reads
