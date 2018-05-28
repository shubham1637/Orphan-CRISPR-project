import random
from random import randint
import re
import subprocess

ConstOffSet = 2000

print('sss')
def findORFStartStop(FileName):
    ORFRanges = []
    for line in open(FileName, "r"):
        RangeMatch = re.search(r'range\s([0-9]{1,})..([0-9]{1,})', line)
        if RangeMatch:
           # LineValues = Line[:-1].split(" ")
            range = str(RangeMatch.group(0))[6:]
            range = [int(i) for i in range.split("..")]
            if range[0] < range[1]:
                ORFRanges.append([range[0] + int(RandomSeqStart) - ConstOffSet, range[1] + int(RandomSeqStart) - ConstOffSet])  # но для обратных рамок будет SeqEnd - range[0] (!) - ConstOFFSet, SeqEnd-range[1] - ConstOFF
            else:
                ORFRanges.append([int(RandomSeqStart + int(LineValues[5])) - range[0] + ConstOffSet, int(RandomSeqStart + int(LineValues[5])) - range[1] + ConstOffSet])
    return ORFRanges

def calculateCoverage(CrisprSeqStart, CrisprSeqEnd, ORFSeqStart, ORFSeqEnd):
    CoverageType = []
    if ORFSeqStart < ORFSeqEnd:  # frames  1, 2, 3
        if (CrisprSeqStart > ORFSeqEnd) | (CrisprSeqEnd < ORFSeqStart) | (
            (CrisprSeqStart < ORFSeqStart) & (CrisprSeqEnd > ORFSeqEnd)):
            CoverageType.append([0, 'None'])
        if (CrisprSeqStart < ORFSeqStart) & (CrisprSeqEnd < ORFSeqEnd) & (CrisprSeqEnd > ORFSeqStart):
            Coverage = round(float(CrisprSeqEnd - ORFSeqStart) / float(CrisprSeqEnd - CrisprSeqStart)  , 2)
            CoverageType.append([Coverage, 'Partial'])
        if (CrisprSeqStart > ORFSeqStart) & (CrisprSeqEnd > ORFSeqEnd) & (CrisprSeqStart < ORFSeqEnd):
            Coverage = round(float(ORFSeqEnd - CrisprSeqStart) / float(CrisprSeqEnd - CrisprSeqStart), 2)
            CoverageType.append([Coverage, 'Partial'])
        if (CrisprSeqStart > ORFSeqStart) & (CrisprSeqEnd < ORFSeqEnd):
            CoverageType.append([1, 'Full'])
    return CoverageType



ORFResultsFileName = "tmp.txt"
# LeftCRISPRsFileName = "/home/utkinai2/Project1/lefties.csv"
IdentifiedCRISPRsFileName = "/home/utkinai2/Project1/identified.csv"
FileWithIDandContigs = "/home/utkinai2/Project1/all1603.pp.txt"
IdentifiedRandomsWithORFFileName = "/home/utkinai2/Project1/Randoms_with_ORF.txt"
RandomSeqCrisprLengthFileName = "/home/utkinai2/Project1/RandomSeq_with_CRISPRlength.txt"
# FileName = 'C:/Users/utkinai2/Desktop/Ipynb-scripts/lefties.csv'

ContigIDSizes = {}
for line in open(FileWithIDandContigs, "r"):
    LineValues = line[:-1].split("\t")
    if line[0] == "#":
        continue
    ContigIDSizes[LineValues[1]] = int(LineValues[4])

count = 0
with open(IdentifiedRandomsWithORFFileName, "w") as IdentifiedCRISPRsWithORFFile:
    for Line in open(IdentifiedCRISPRsFileName, "r"):
        count += 1
        if count < 2:
            continue         # header is skipped now
        LineValues = Line[:-1].split(",")  # the last one is perenos stroki, thats'why -1

        RandomSeqStart = randint(ConstOffSet+1, max(ConstOffSet +2,ContigIDSizes[LineValues[2].strip("\"")]-ConstOffSet-int(LineValues[5]) -1))
        subprocess.call("blastdbcmd -db /panfs/pan1/prokdata/db/all1603.nt" + " -entry " + LineValues[2].strip("\"") +
                        " -range " + str(RandomSeqStart - ConstOffSet) + "-" + str(RandomSeqStart + int(LineValues[5])+ ConstOffSet) + " | orf -f=0 -m=1 -n=1 > " + ORFResultsFileName , shell=True)
#     # for Line in open(ORFResultsFileName, "r"):
#     #     LeftCRISPRsWithORFFile.write(Line + "\n")

    #ORFranges = findORFStartStop(LeftCRISPRsWithORFFile)
        ORFranges = findORFStartStop(ORFResultsFileName)
        # print(ORFranges)
        CoverageTypesArrays = []
        for range in ORFranges:
            for elem in calculateCoverage(int(RandomSeqStart), int(RandomSeqStart + int(LineValues[5])), range[0], range[1]):
                CoverageTypesArrays.append([elem, range[0], range[1]])
            #print(CoverageTypesArrays)

        BestCoverage = 0
        BestORFType = []
        for element in CoverageTypesArrays:
            if BestCoverage <= element[0][0]:
                BestCoverage = element[0][0]
                BestORFType = element
        print(BestORFType)
        #print(BestCoverage)
        IdentifiedCRISPRsWithORFFile.write(LineValues[2].strip("\"") + ' ' + str(BestORFType[0][0]) + ' ' + BestORFType[0][1] + ' ' + str(BestORFType[1]) + ' ' + str(BestORFType[2]) + ' ' +
                                     str(RandomSeqStart)+ ' ' + str(RandomSeqStart + int(LineValues[5])) + "\n")



