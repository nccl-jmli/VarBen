from subprocess import call
import pysam
from varben.common.methods import getReadStrand
import random
import os


def _call(cmd):
    print cmd
    flag = call(cmd, shell=True)
    if flag == 0:
        return True
    else:
        raise Exception("Cmd Error: %s" % cmd)


def bamToFastq(bamFile, outPrefix, is_single):
    if not is_single:
        fq1 = outPrefix + "_1.fq"
        fq2 = outPrefix + "_2.fq"
        bamToFastq_cmd = "bedtools bamtofastq  -i %s -fq %s -fq2 %s" % (bamFile, fq1, fq2)
    else:
        fq1 = outPrefix + ".fq"
        fq2 = None
        bamToFastq_cmd = "bedtools bamtofastq  -i %s -fq %s " % (bamFile, fq1)
    print "bam to fastq start ..................................."
    _call(bamToFastq_cmd)
    print "bam to fastq end ....................................."
    return fq1, fq2


def bamToFastq_2(bamFile, outPrefix, is_single):
    fq1 = outPrefix + ".fq"
    bamToFastq_cmd = "bedtools bamtofastq  -i %s -fq %s " % (bamFile, fq1)
    print "bam to fastq start ..................................."
    _call(bamToFastq_cmd)
    print "bam to fastq end ....................................."
    return fq1


def map_bwa(ref_index, outSamFile, fq1, fq2=None, header=None, threadNum=1):
    RG_head = ''
    if header:
        RG_head = '-R "%s"' % header
    if fq2 is None:
        mapping_cmd = "bwa mem %s -t %s %s %s >%s" % (RG_head, threadNum, ref_index, fq1, outSamFile)
    else:
        mapping_cmd = "bwa mem %s -t %s %s %s %s >%s" % (RG_head, threadNum, ref_index, fq1, fq2, outSamFile)
    print "mapping by bwa start ...................................."
    _call(mapping_cmd)
    print "mapping by bwa end ......................................"


def map_novoalign(ref_index, out_sam, fq1, fq2=None, process=1):
    if fq2 is None:
        # mapping_cmd = "/lustre/users/fangshuangsang/Project/simulation_mutation/softs/novocraft/novoalign -d %s -f %s --mCPU %s -o SAM > %s" % (ref_index, fq1, process, out_sam)
        mapping_cmd = "novoalign -d %s -f %s --mCPU %s -o SAM > %s" % (ref_index, fq1, process, out_sam)
    else:
        # mapping_cmd = "/lustre/users/fangshuangsang/Project/simulation_mutation/softs/novocraft/novoalign -d %s -f %s %s --mCPU %s -o SAM > %s" % (ref_index, fq1, fq2, process, out_sam)
        mapping_cmd = "novoalign -d %s -f %s %s --mCPU %s -o SAM > %s" % (ref_index, fq1, fq2, process, out_sam)
    print "mapping by novoalign start ...................................."
    _call(mapping_cmd)
    print "mapping by novoalign end ...................................."


def samToBam(inSamFile, outBamFile):
    samToBam_cmd = "samtools view -bS %s > %s" % (inSamFile, outBamFile)
    _call(samToBam_cmd)


# def bamSort_new(inBamFile, sortedBamFile, sort_key=None):
def bamSort(inBamFile, sortedBamFile_prefix, sort_key=None):
    if not sort_key:
        bamSort_cmd = "samtools sort -o %s %s" % (sortedBamFile_prefix + ".bam", inBamFile)
    elif sort_key == "name":
        bamSort_cmd = "samtools sort -n -o %s %s" % (sortedBamFile_prefix + ".bam", inBamFile)
    print "samtools sort start................................"
    _call(bamSort_cmd)
    print "samtools sort end.................................."


def bamIndex(inBamFile):
    index_cmd = "samtools index %s" % (inBamFile)
    print "samtools index start .............................."
    _call(index_cmd)
    print "samtools index end ................................"


def remap(ref_index, inBamFile, outBamFile, aligner, is_single, header=None, sort=True, threadNum=4):
    print "remap start ......................................."
    aligner = aligner.lower()
    prefix = outBamFile.rstrip(".bam")
    if aligner in ("bwa", "novoalign"):
        # prefix = outBamFile.rstrip(".bam")
        if sort:
            inPrefix = inBamFile.rstrip(".bam")
            bam_toConvert = inPrefix + ".sortByName.bam"
            bamSort(inBamFile, inPrefix + ".sortByName", sort_key="name")
        else:
            bam_toConvert = inBamFile
        fastqPrefix = prefix + ".to"
        fq1, fq2 = bamToFastq(bam_toConvert, fastqPrefix, is_single)
        outSamFile = prefix + ".sam"

        if aligner == "bwa":
            map_bwa(ref_index, outSamFile, fq1, fq2, header=header, threadNum=threadNum)
        elif aligner == "novoalign":
            map_novoalign(ref_index, outSamFile, fq1, fq2, threadNum)
        outBamFile_tmp = prefix + ".unsorted.bam"
        samToBam(outSamFile, outBamFile_tmp)
        bamSort(outBamFile_tmp, prefix)
        bamIndex(outBamFile)
    elif aligner == "tmap":
        outBamFile_tmp = prefix + ".unsorted.bam"
        remap_tmap(ref_index, inBamFile, outBamFile_tmp, threadNum)
        bamSort(outBamFile_tmp, prefix)
        bamIndex(outBamFile)
    print "remap end .........................................."
    return outBamFile


# discard
def remap_2(ref_index, inBamFile, outBamFile, aligner, is_single, threadNum=1):
    print "remap start ......................................."
    if aligner == "bwa":
        prefix = outBamFile.rstrip(".bam")
        fastqPrefix = prefix + ".to"
        fq = bamToFastq_2(inBamFile, fastqPrefix, is_single)
        outSamFile = prefix + ".sam"
        map_bwa(ref_index, outSamFile, fq, threadNum=threadNum)
        outBamFile_tmp = prefix + ".unsorted.bam"
        samToBam(outSamFile, outBamFile_tmp)
        bamSort(outBamFile_tmp, prefix)
        bamIndex(outBamFile)
    print "remap end .........................................."
    return outBamFile


def remap_tmap(ref_index, inBamFile, outBamFile, threadNum=4):
    mapping_cmd = "tmap mapall -n %s -f %s -r %s -i bam -s %s -o 2 -v -Y -u --prefix-exclude 5 -o 2 -J 25 --context stage1 map4" % (
        threadNum, ref_index, inBamFile, outBamFile)
    _call(mapping_cmd)
    return outBamFile


def bamMerge(bamList, outBamFile):
    cmd = "samtools merge -c -f %s %s" % (outBamFile, " ".join(bamList))
    _call(cmd)


def bamAddRG_picard(inBamFile, outBamFile, picard_path, RGID=1, RGLB='BE', RGPL='illumina', RGPU='BE', RGSM='BE'):
    cmd = "java -jar %s/AddOrReplaceReadGroups.jar I=%s O=%s RGID=%s RGLB=%s \
    RGPL=%s RGPU=%s RGSM=%s SO=coordinate TMP_DIR=%s" % (
        picard_path, inBamFile, outBamFile, RGID, RGLB, RGPL, RGPU, RGSM,
        os.path.join(os.path.dirname(outBamFile), 'tmp'))
    _call(cmd)


def bamMerge_picard(bamList, outBamFile, picard_path, sort_order='coordinate'):
    cmd = "java -jar %s/MergeSamFiles.jar SORT_ORDER=%s %s OUTPUT=%s TMP_DIR=%s & " % (
        picard_path, sort_order, " ".join(["INPUT=" + bamfile for bamfile in bamList]), outBamFile,
        os.path.join(os.path.dirname(outBamFile), 'tmp'))
    _call(cmd)


def getRegionReads(inBam, regionBed, inRegionBam, outRegionBam):
    cmd = "samtools view %s -b -h -o %s -U %s -L %s" % (inBam, inRegionBam, outRegionBam, regionBed)
    _call(cmd)


def bamAddRG(editRemap, editBamReads, templateBamFile, outBamFile):
    # editRemapBam_addRG_File = tempOutDir + "/edit.remap.addRG.bam"
    head = editRemap.header
    head["RG"] = templateBamFile.header["RG"]
    addRGBam = pysam.AlignmentFile(outBamFile, 'wb', header=head)
    RG = _getRGs(templateBamFile)
    for read in editRemap.fetch():
        readName = read.query_name
        strand = getReadStrand(read)
        if readName in editBamReads:
            orig = editBamReads[readName][strand]
        else:
            orig = None
        newRead = readAddRG(read, orig, RG)
        # print newRead
        addRGBam.write(newRead)
    addRGBam.close()


def readAddRG(read, orig, RG=None):
    if RG:
        hasRG = False
        if read.tags is not None:
            for tag in read.tags:
                if tag[0] == 'RG':
                    hasRG = True

        # use RG from original read if it exists
        if orig is not None:
            if not hasRG and orig.tags is not None:
                for tag in orig.tags:
                    if tag[0] == 'RG':
                        read.tags = read.tags + [tag]
                        hasRG = True

        if not hasRG:
            # give up and add random read group from list in header (e.g. for simulated reads)
            newRG = RG[random.randint(0, len(RG) - 1)]
            read.tags = read.tags + [("RG", newRG)]
    return read


def _getRGs(bam):
    '''return list of RG IDs'''
    RG = []
    if 'RG' in bam.header:
        for headRG in bam.header['RG']:
            RG.append(headRG['ID'])
    return RG


def bamReadAddTag(bamFile, tag, outFile):
    bam = pysam.AlignmentFile(bamFile, 'rb')
    outBam = pysam.AlignmentFile(outFile, 'wb', template=bam)
    for read in bam.fetch():
        tags = read.tags
        tags.append((tag, 1))
        read.tags = tags
        outBam.write(read)
    outBam.close()
    bam.close()
