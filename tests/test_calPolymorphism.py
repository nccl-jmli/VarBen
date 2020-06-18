import pysam

bam = pysam.AlignmentFile('./U2AF1_sorted.bam',"rb")
chrom='chr21'
start=44524455-1
end=44524455
reads = bam.fetch(chrom, start, end + 1)
numDict = {}
readsNum = 0
for pos in range(start, end):
    numDict[pos] = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    
for read in bam.fetch(chrom, start, end + 1):
    readsNum += 1
    #print(read)
    for pair in read.get_aligned_pairs():
        if start <= pair[1] < end:
            pos = pair[0]
            if pos is None:
                continue
            if pos < len(read.query_alignment_sequence): # should be 'query_length'
                nul = read.query_alignment_sequence[pos] # should be 'query_sequence'
                #if nul in 'ATCG':
                #    numDict[pair[1]][nul] += 1
                if nul in 'TCG' and pos < 10:
                    print( "[POSITON: "+str(pos) +"]\t [BASE: "+ nul+ "]\t [READ: "+str(read)+"]\t READPOS: "+str(read.get_aligned_pairs())+"]\t READBASE:"+str(read.query_alignment_sequence[pos])+"]")
                    break
    if nul in 'TCG' and pos < 10:
        break
