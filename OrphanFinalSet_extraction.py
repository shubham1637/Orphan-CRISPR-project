import re

OrphansFileName = "/home/utkinai2/Project1/Orphans.txt"
OrphanSetFileName = "/home/utkinai2/Project1/OrphanSet_with_uniqueSpacers.txt"
OrphanWithUniqueSpacersFileName = "/home/utkinai2/Project1/Orphans_with_uniqueSpacers.txt"
OrphansWithProtHitsFileName = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/orphancrispr/IDorphans_withHits_inProt.csv"
ProteinSetFileName = "/home/utkinai2/Project1/Orphans_withProtHits.txt"
HitsFileName = "/home/utkinai2/Project1/AllOrphanHits_inProt.txt"
FinalFileName = "/home/utkinai2/Project1/Orphans_finalfiltered.txt"

Orphans = {}
# with open(OrphanWithUniqueSpacersFileName, "w") as f:
for line in open(OrphanSetFileName, "r"):
    LineValues = line[:-1].split(" ")
    if LineValues[1].replace('"', '') not in Orphans:
        Orphans[LineValues[1].replace('"', '')] = [[LineValues[2], LineValues[3]]]
    else:
        Orphans[LineValues[1].replace('"', '')] += [[LineValues[2], LineValues[3]]]

    # for line in open(OrphansFileName,"r"):
    #     LineValues = line[:-1].split(" ")
    #     if LineValues[2].replace('"', '') in Orphans:
    #         for case in Orphans[LineValues[2].replace('"', '')]:
    #             if case == [LineValues[3], LineValues[4]]:
    #                 f.write(line)

count = 0
# for line in open(OrphanWithUniqueSpacersFileName, "r"):
for line in open(OrphanWithUniqueSpacersFileName, "r"):
    count += 1
print(count)
count3 = 0
# getting rid of orphan with meaningful hits in proteins
BadSet = {}
with open(HitsFileName, "w") as HitsFile:
    with open(ProteinSetFileName, "w") as ProtFile:
        for line in open(OrphansWithProtHitsFileName, "r"):
            LineValues = line[:-1].split(",")
            ID = LineValues[0].split("_")
            if ID[0] in Orphans:
                for case in Orphans[ID[0]]:
                    if case == [ID[1], ID[2]]:
                        if ID[0] not in BadSet:
                            BadSet[ID[0]] = [[ID[1], ID[2], LineValues[-3]]]
                        else:
                            if ID[1] !=  BadSet[ID[0]][0][0]:
                                BadSet[ID[0]] += [[ID[1], ID[2],  LineValues[-3]]]
            HitsFile.write(ID[0] + '\t' + str(ID[1]) + '\t' + str(ID[2]) + '\t' + LineValues[-3] + '\n')
        for key in BadSet:
            for badcase in BadSet[key]:
                count3 += 1
                ProtFile.write(key + '\t' + str(badcase[0]) + '\t' + str(badcase[1]) + '\t' + badcase[2] + '\n')
print(count3)

count2 = 0
with open(FinalFileName, "w") as GoodFile:
    for line in open(OrphanWithUniqueSpacersFileName, "r"):
        LineValues = line[:-1].split(" ")
        if LineValues[2].replace('"', '') in BadSet:
            flag = count2
            for case in BadSet[LineValues[2].replace('"', '')]:
                if [case[0], case[1]] == [LineValues[3], LineValues[4]]:
                    count2 += 1
            if flag == count2:
                GoodFile.write(line)
        else:
            GoodFile.write(line)
print(count2)

count1 = 0
for line in open(FinalFileName, "r"):
    count1 += 1
print(count1)