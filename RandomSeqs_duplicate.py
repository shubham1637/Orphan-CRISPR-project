import random
from random import randint
import re
import subprocess

ConstOffSet = 2000

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
                ORFRanges.append([int(RandomSeqStart + int(LineValues[6])) - range[0] + ConstOffSet, int(RandomSeqStart + int(LineValues[6])) - range[1] + ConstOffSet])
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
FileWithIDandContigs = "/panfs/pan1/prokdata/db_tmp/all1603.pp.txt"
RandomsWithORFFileName = "/home/utkinai2/Project1/Randoms_with_ORF_full-part.txt"

# FileName = 'C:/Users/utkinai2/Desktop/Ipynb-scripts/lefties.csv'

ContigIDSizes = {}
for line in open(FileWithIDandContigs, "r"):
    LineValues = line[:-1].split("\t")
    if line[0] == "#":
        continue
    ContigIDSizes[LineValues[1]] = int(LineValues[4])

count = 0
with open(RandomsWithORFFileName, "w") as RandomsWithORFFile:
    for Line in open(IdentifiedCRISPRsFileName, "r"):
        count += 1
        if count < 2:
            continue         # header is skipped now
        LineValues = Line[:-1].split(",")  # the last one is perenos stroki, thats'why -1
        Rand1 = randint(1, max(1, int(LineValues[3])-int(LineValues[5])))
        Rand2 = randint(min(int(LineValues[4]), int(ContigIDSizes[LineValues[2].strip("\"")]) - int(LineValues[5])),int(ContigIDSizes[LineValues[2].strip("\"")]) - int(LineValues[5]))
        RandomSeqStart = random.choice([Rand1,Rand2])
        subprocess.call("blastdbcmd -db /panfs/pan1/prokdata/db/all1603.nt" + " -entry " + LineValues[2].strip("\"") +
                        " -range " + str(RandomSeqStart) + "-" + str(RandomSeqStart + int(LineValues[5])) + " | orf -f=0 -m=1 -n=1 > " + ORFResultsFileName , shell=True)

        # # if orf doesn't work - export PATH=$PATH:~wolf/bin/ and then do: export PATH=$PATH:~wolf/bin/add/

        MaxORFlength = 0
        for line in open(ORFResultsFileName, "r"):
            if line[0] == ">":
                continue
            MaxORFlength = max(MaxORFlength, len(line[:-1]) * 3)
        if MaxORFlength < int(LineValues[4]) - int(LineValues[3]) - 3:
            RandomsWithORFFile.write(
                LineValues[2] + ' ' + str(LineValues[3]) + ' ' + str(LineValues[4]) + ' ' + "Partial" + ' ' + str(
                    int(LineValues[4]) - int(LineValues[3])) + '\n')
        else:
            RandomsWithORFFile.write(
                LineValues[2] + ' ' + str(LineValues[3]) + ' ' + str(LineValues[4]) + ' ' + "Full" + ' ' + str(
                    int(LineValues[4]) - int(LineValues[3])) + '\n')
        print(count)
    #ORFranges = findORFStartStop(LeftCRISPRsWithORFFile)
        # ORFranges = findORFStartStop(ORFResultsFileName)
        # # print(ORFranges)
        # CoverageTypesArrays = []
        # for range in ORFranges:
        #     for elem in calculateCoverage(int(RandomSeqStart), int(RandomSeqStart + int(LineValues[6])), range[0], range[1]):
        #         CoverageTypesArrays.append([elem, range[0], range[1]])
        #     #print(CoverageTypesArrays)
        #
        # BestCoverage = 0
        # BestORFType = []
        # for element in CoverageTypesArrays:
        #     if BestCoverage <= element[0][0]:
        #         BestCoverage = element[0][0]
        #         BestORFType = element
        # print(BestORFType)
        # #print(BestCoverage)
        # BigIdentifiedCRISPRsWithORFFile.write(LineValues[3].strip("\"") + ' ' + str(BestORFType[0][0]) + ' ' + BestORFType[0][1] + ' ' + str(BestORFType[1]) + ' ' + str(BestORFType[2]) + ' ' +
        #                              str(RandomSeqStart)+ ' ' + str(RandomSeqStart + int(LineValues[6])) + "\n")
