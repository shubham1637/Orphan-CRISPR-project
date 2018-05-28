import re
from collections import defaultdict
import subprocess

# FileName = '/home/utkinai2/Project1/lefties.csv'
GoodRange = defaultdict(list)

ConstOffSet = 2000
LineValuesORF = ()
# FileNameORF = '/home/utkinai2/OutputTest1.txt'
LeftCRISPRsWithORFFile = 'C:/Users/utkinai2/Desktop/Ipynb-scripts/OutputTest134.txt'
FileName = 'C:/Users/utkinai2/Desktop/Ipynb-scripts/lefties.csv'
count = 0


def findORFStartStop(File):
    ORFRanges = []
    Lines = [line[:-1].rstrip(' ') for line in open(File, "r")]
    for i in Lines:
        LineValuesORF = str(i)
        RangeMatch = re.search(r'range\s([0-9]{1,})..([0-9]{1,})', LineValuesORF)
        if RangeMatch:
            range = str(RangeMatch.group(0))[6:]
            range = [int(i) for i in range.split("..")]
            ORFRanges.append([range[0] + SeqStart, range[1] + SeqStart])
    return ORFRanges

def calculateCoverage(CrisprSeqStart, CrisprSeqEnd, ORFSeqStart, ORFSeqEnd):
    CoverageType = []
    if ORFSeqStart < ORFSeqEnd:  # frames  1, 2, 3
        if (CrisprSeqStart > ORFSeqEnd) | (CrisprSeqEnd < ORFSeqStart) | (
            (CrisprSeqStart < ORFSeqStart) & (CrisprSeqEnd > ORFSeqEnd)):
            CoverageType.append([0, 'None'])
        if (CrisprSeqStart <= ORFSeqStart) & (CrisprSeqEnd <= ORFSeqEnd) & (CrisprSeqEnd >= ORFSeqStart):
            Coverage = round((CrisprSeqEnd - ORFSeqStart) / (CrisprSeqEnd - CrisprSeqStart), 1)
            print(Coverage)
            CoverageType.append([Coverage, 'right'])
        if (CrisprSeqStart > ORFSeqStart) & (CrisprSeqEnd > ORFSeqEnd) & (CrisprSeqStart < ORFSeqEnd):
            Coverage = round((ORFSeqEnd - CrisprSeqStart) / (CrisprSeqEnd - CrisprSeqStart), 1)
            CoverageType.append([Coverage, 'left'])
        if (CrisprSeqStart > ORFSeqStart) & (CrisprSeqEnd < ORFSeqEnd):
            Coverage = round((ORFSeqEnd - ORFSeqStart) / (CrisprSeqEnd - CrisprSeqStart), 1)
            CoverageType.append([Coverage, 'inside'])
    if ORFSeqStart > ORFSeqEnd:  # frames  -1, -2, -3
        if (CrisprSeqStart > ORFSeqStart) | (CrisprSeqEnd < ORFSeqEnd) | (
            (CrisprSeqStart < ORFSeqEnd) & (CrisprSeqEnd > ORFSeqStart)):
            CoverageType.append([0, 'None'])
        if (CrisprSeqStart < ORFSeqEnd) & (CrisprSeqEnd < ORFSeqStart) & (CrisprSeqEnd > ORFSeqEnd):
            Coverage = round((CrisprSeqEnd - ORFSeqEnd) / (CrisprSeqEnd - CrisprSeqStart), 1)
            CoverageType.append([Coverage, 'right'])
        if (CrisprSeqStart > ORFSeqEnd) & (CrisprSeqEnd > ORFSeqStart) & (CrisprSeqStart < ORFSeqStart):
            Coverage = round((ORFSeqStart - CrisprSeqStart) / (CrisprSeqEnd - CrisprSeqStart), 1)
            CoverageType.append([Coverage, 'left'])
        if (CrisprSeqStart > ORFSeqEnd) & (CrisprSeqEnd < ORFSeqStart):
            Coverage = round((ORFSeqStart - ORFSeqEnd) / (CrisprSeqEnd - CrisprSeqStart), 1)
            CoverageType.append([Coverage, 'inside'])
    return CoverageType

for Line in open(FileName, "r"):
    count += 1
    if count < 135:
        continue     # header is skipped now
    LineValues = Line[:-1].split(",")  # the last one is perenos stroki, thats'why -1
    SeqStart = int(LineValues[3])  # Вручную вычтены 2000 (и прибавлены для следующей строки) при получении текстового файла с ОРФами для теста
    SeqEnd = int(LineValues[4])
    print(SeqStart, ' ' ,SeqEnd)

    #print(SeqEnd , ' ',  SeqStart)
    # но в текстовом файле (ниже) будет диапазон -2000 + SeqStart --- SeqEnd +2000 !!!!!!!!!!!! KEEP IN MIND

    # subprocess.call("blastdbcmd -db /panfs/pan1/prokdata/db/all1603.nt" + " -entry " + LineValues[2] + " -range " +
    #                 SeqStart + "-" + SeqEnd, shell=True)


# This part checks the length of ORF, then if length of array is shorter, it saves its ID and checks
# whether the found ORF covers partially or fully the given array (through range parsing)

    ORFranges = findORFStartStop(LeftCRISPRsWithORFFile)
    print(ORFranges)

    for range in ORFranges:
        print(calculateCoverage(SeqStart, SeqEnd, range[0], range[1]))
        print(SeqStart,' ', SeqEnd,' ', range[0],' ', range[1])
    break

    #if LengthMatch:
        #length = LengthMatch.group(1)
        # SeqMatch = re.findall(r'(CP[0-9]{5,}.[0-9])', LineValuesORF)
        # for i in SeqMatch:
        #SeqBoundaries = [SeqStart + range[0] + ConstOffSet, SeqStart + range[1] + ConstOffSet]
        # if abs(SeqBoundaries[1] - SeqBoundaries[0]) > MaxORF:
        #print(dict(GoodRange))
#             if (range[1] > range[0])&(SeqStart < SeqStart - 2000 + range[0] < SeqEnd)&(SeqStart - 2000 + range[1] > SeqEnd):
#             # frames 1, 2, 3, partially covered, ORF covers only SeqEnd
#                 GoodRange['right'].append(SeqBoundaries)
#             if (range[1] > range[0]) & (range[0] < 2000) & (SeqStart - 2000 + range[1] < SeqEnd):
#             # frames 1, 2, 3, partially covered, ORF covers only SeqStart
#                 GoodRange['left'].append(SeqBoundaries)
#             if (range[1] > range[0])&(range[0] < 2000)&(SeqStart - 2000 + range[1] > SeqEnd):
#             # frames 1, 2, 3, fully covered
#                 GoodRange['inside'].append(SeqBoundaries)
#             if (range[1] < range[0])&(SeqStart < SeqStart - 2000 + range[1] < SeqEnd)&(SeqStart - 2000 + range[0] > SeqEnd):
#             # frames -1, -2, -3, partially covered, ORF starts between SeqStart and SeqEnd
#                 GoodRange['right'].append(SeqBoundaries)
#             if (range[1] < range[0]) & (range[1] < 2000) & (SeqStart - 2000 + range[0] < SeqEnd):
#             # frames -1, -2, -3, partially covered, ORF ends between SeqStart and SeqEnd
#                 GoodRange['left'].append(SeqBoundaries)
#             if (range[1] < range[0])&(range[1] < 2000)&(SeqStart - 2000 + range[0] > SeqEnd):
#             # frames -1, -2, -3, fully covered
#                 GoodRange['inside'].append(SeqBoundaries)
# print(list(GoodORFs.items())[-1][1])
# BestBoundaries = list(GoodORFs.items())[-1][1]


# ListOfMappedORFs = str()
# ListOfGoodORFs = str()
# for i in list(GoodORFs.values()):
#     ListOfGoodORFs += str(i)+ ","
# print(ListOfGoodORFs)
# for i in GoodRange.values():
#     for j in i:
#         ListOfMappedORFs += str(j) + ","
# print(ListOfMappedORFs)

# subprocess.call(
#     "blastdbcmd -db " + Database + " -entry_batch " + ClusterGIListFileName + " > " + ClusterFASTA, shell=True)
#