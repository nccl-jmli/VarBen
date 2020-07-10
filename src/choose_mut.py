#!/usr/bin/env python
import sys
import random

def choose(fileName, chooseNum, outFile):
    lineList = []
    with open(fileName, 'r') as fin:
        lineList = fin.readlines()

    fin.close()
    print len(lineList), type(lineList)
    randomLines = random.sample(lineList, chooseNum)
    with open(outFile, 'w') as fout:
        for line in randomLines:
            fout.write(line)
    fout.close()
    
if len(sys.argv) == 4:
    fileName = sys.argv[1]
    chooseNum = int(sys.argv[2])
    outFile = sys.argv[3]
    print fileName, chooseNum, outFile
    choose(fileName, chooseNum, outFile)
else:
    sys.exit("usage: %s <filename> <Choose number> <output filename with path> " % sys.argv[0])

