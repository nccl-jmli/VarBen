__author__ = 'fangshuangsang'


class Read:
    def __init__(self, read):
        self.query_name = read.query_name
        self.strand = self.getReadStrand(read)
        self.keyName = self.getKeyName()
        self.reference_name = read.reference_name
        self.query_sequence = read.query_sequence
        self.query_qualities = read.query_qualities
        self.quality = "".join([chr(x + 33) for x in self.query_qualities])
        self.get_aligned_pairs = read.get_aligned_pairs()
        self.reference_start = read.reference_start
        self.reference_end = read.reference_end
        self.is_reverse = read.is_reverse
        self.is_read1 = read.is_read1
        self.is_read2 = read.is_read2

    def __str__(self):
        return self.keyName + "\t" + self.reference_name + "\t" + self.query_sequence + "\t" + str(self.query_qualities) + "\n"

    def getReadStrand(self, read):
        if read.is_read1:
            return '1'
        elif read.is_read2:
            return '2'
        else:
            return '0'

    def getKeyName(self):
        return self.query_name + "/" + self.strand


class Mutation:
    def __init__(self, chrom, start, end, freq, muttype, alt):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.freq = freq
        self.muttype = muttype
        self.alt = alt

    def __str__(self):
        return "chrom: %s, start: %s, end: %s, freq: %s, muttype: %s, alt: %s" % (
            self.chrom, self.start, self.end, self.freq, self.muttype, self.alt)

    def mutinfo(self):
        return "%s\t%s\t%s\t%s\t%s\t%s" % (
            self.chrom, self.start + 1, self.end + 1, self.freq, self.muttype, self.alt)


class Haplotype(object):
    def __init__(self, haplotype_chrom, haplotype_start, haplotype_end, haplotype_freq, mutList):
        self.chrom = haplotype_chrom
        self.start = haplotype_start
        self.end = haplotype_end
        self.freq = haplotype_freq
        self.mutList = mutList

    def __str__(self):
        return "chrom: %s, start: %s, end: %s, freq: %s, mutList:[ %s ]" % (
            self.chrom, self.start, self.end, self.freq, ";".join([str(mut) for mut in self.mutList]))

    def mutinfo(self):
        return "\n".join([mut.mutinfo() for mut in self.mutList])


class SV(object):
    def __init__(self, chrom, start, end, sv_type, freq, trans_chrom=None, trans_start=None, trans_end=None, cnv_type=None, dup_num=None):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.freq = float(freq)
        self.sv_type = sv_type
        self.trans_chrom = trans_chrom
        self.trans_start = trans_start
        self.trans_end = trans_end
        self.cnv_type = cnv_type
        self.dup_num = dup_num

    def __str__(self):
        str_sv = "chrom: %s, start: %s, end: %s, sv_type: %s, freq: %s" % (
            self.chrom, self.start, self.end, self.sv_type, self.freq)
        if self.trans_chrom:
            str_sv = str_sv + ", trans_chrom: %s, trans_start: %s, trans_end: %s" % (
                self.trans_chrom, self.trans_start, self.trans_end)
        elif self.cnv_type:
            str_sv = str_sv + ", cnv_type: %s" % self.cnv_type
        elif self.dup_num:
            str_sv = str_sv + ", dup_num: %s" % self.dup_num
        return str_sv

    def svinfo(self):
        str_sv = "%s\t%s\t%s\t%s\t%s" % (
            self.chrom, self.start + 1, self.end + 1, self.sv_type, self.freq)
        if self.trans_chrom:
            str_sv = str_sv + "\t%s\t%s\t%s" % (
                self.trans_chrom, self.trans_start + 1, self.trans_end + 1)
        elif self.cnv_type:
            str_sv = str_sv + "\t%s" % self.cnv_type
        elif self.dup_num:
            str_sv = str_sv + "\t%s" % self.dup_num
        return str_sv


