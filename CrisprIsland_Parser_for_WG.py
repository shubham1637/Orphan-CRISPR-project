import re

BothInWGFileName = "/home/utkinai2/Project1/AllOrphans_and_Id_inWG.txt"
IslandsAnnotationFileName = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/CRISPR_islands_new/Islands.ann_CRISPR"
WGFileName = "/home/utkinai2/Project1/near_orphans_inWG.txt"

BothIDandRange = {}
OrphansIDandRange = {}
count = 0
Const = 2000

for line in open(BothInWGFileName, "r"):
    LineValues = line[:-1].split(" ")
    IDMatch = re.search(r'([A-Z]{2,}([0-9]{1,}))', LineValues[2])
    if IDMatch:
        if str(IDMatch.group(0)) not in BothIDandRange:
            BothIDandRange[str(IDMatch.group(0))] = [[LineValues[3],LineValues[4], LineValues[1].replace('"', ''), LineValues[6].replace('"', '')]]
        else:
            BothIDandRange[str(IDMatch.group(0))] += [[LineValues[3], LineValues[4], LineValues[1].replace('"', ''),
                                                     LineValues[6].replace('"', '')]]
# print(BothIDandRange)

flag = 0
with open(WGFileName, "w") as WGFile:
    for line in open(IslandsAnnotationFileName, "r"):
        print(flag)
        LineValues = line[:-1].split("\t")  # в строках-заголовках [1] - это пробел!!! а в строках обычных [7] - это пробел!!! (wtf, but keep in mind)
        if LineValues[0] == "===":
            ID = LineValues[4]
            if ID in BothIDandRange:
                flag = 1
            else:
                flag = 0
        else:
            if flag == 1:
                for i in range(len(BothIDandRange[LineValues[4]])):
                    if (int(BothIDandRange[LineValues[4]][i][0]) > (int(re.search(r'([0-9]{1,})..([0-9]{1,})', LineValues[1]).group(0).split("..")[0]) - 10000)) and \
                        (int(BothIDandRange[LineValues[4]][i][0]) < ((int(re.search(r'([0-9]{1,})..([0-9]{1,})', LineValues[1]).group(0).split("..")[1])) + 10000)) or \
                        (int(BothIDandRange[LineValues[4]][i][1]) > (int(re.search(r'([0-9]{1,})..([0-9]{1,})', LineValues[1]).group(0).split("..")[0]) - 10000)) and \
                             (int(BothIDandRange[LineValues[4]][i][1]) < (int(re.search(r'([0-9]{1,})..([0-9]{1,})', LineValues[1]).group(0).split("..")[1]) + 10000)):
                        if LineValues[8] != "Unknown":
                            WGFile.write(str(LineValues[4]) + '\t' + str(BothIDandRange[LineValues[4]][i][0]) + '\t' + str(BothIDandRange[LineValues[4]][i][1]) + '\t' +
                            str(BothIDandRange[LineValues[4]][i][2]) + '\t'+ str(BothIDandRange[LineValues[4]][i][3]) + '\t' + str(LineValues[8]) + '\n')

