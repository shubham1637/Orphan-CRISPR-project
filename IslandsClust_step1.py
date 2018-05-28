import re


# IslandsAnnotationFileName = "/panfs/pan1/orphancrispr/IslandsCluster/Islands_pfam.ann_clust"
IslandsAnnotationFileName = "/panfs/pan1/prokdata/CRISPRicity/JointClusters/MergedIslands.tsv"
# IslandIDwithCRISPRNeighborsFileName = "/panfs/pan1/orphancrispr/IslandsCluster/Islands_pfam_3NeighborsSet.txt"
IslandIDwithCRISPRNeighborsFileName = "/panfs/pan1/orphancrispr/IslandsCluster/Islands_Identified_3NeighborsSet.txt"
IslandIDwithCRISPRNeighborsFixedFileName = "/panfs/pan1/orphancrispr/IslandsCluster/Islands_Identified_3NeighborsSet_fixedID.txt"
IdentifiedArraysInfoFileName = "/home/utkinai2/Project1/CRISPR_info_known_021718.tsv"


def Unlist(MyList):
    flattened = str()
    for i in MyList:
        if re.search(r',', i):
            flattened += ' ' + i.split(',')[0] + ' ' + i.split(',')[1]
        else:
            flattened += ' ' + i.strip("'")
    return flattened



IslandIDwithGenes = {}
#with open(IslandIDwithUniqueSetLengthAnd6genesFileName, "w") as IslandsInfoFile:
for line in open(IslandsAnnotationFileName, "r"):
    LineValues = line[:-1].split("\t")  # в строках-заголовках [1] - это пробел!!! а в строках обычных [7] - это пробел!!! (wtf, but keep in mind)
    if LineValues[0] == "===":
        GeneList = []
        count = 0
        flag = 0
        UniqueSet = set()
        CRISPRcoordinates = []
        IslandID = ''
    else:
        count += 1
        if LineValues[6] != "Unknown" and LineValues[6] != "CRISPR" :
            GeneList.append(LineValues[6])
            if (flag == 1 and len(UniqueSet) < min(6, 2 + CRISPRcoordinates[0])) \
                    or (flag == 2 and len(UniqueSet) < min(12, min(6, 2 + CRISPRcoordinates[0]) +
                                                                min(6, CRISPRcoordinates[1] - CRISPRcoordinates[0] - 1))) \
                    or ((flag == 3 and len(UniqueSet) < min(18, min(6, 2 + CRISPRcoordinates[0]) +
                                                                min(6, CRISPRcoordinates[1] - CRISPRcoordinates[0] - 1) +
                                                                min(6, CRISPRcoordinates[2] - CRISPRcoordinates[1] - 1 )))) \
                    or (flag == 4 and len(UniqueSet) < min(24,  min(6, 2 + CRISPRcoordinates[0]) +
                                                                min(6, CRISPRcoordinates[1] - CRISPRcoordinates[0] - 1) +
                                                                min(6, CRISPRcoordinates[2] - CRISPRcoordinates[1] - 1 ) +
                                                                min(6, CRISPRcoordinates[3] - CRISPRcoordinates[2] - 1))):
                UniqueSet.add(LineValues[6])
            if IslandID in IslandIDwithGenes:
                IslandIDwithGenes[IslandID] = [list(UniqueSet), len(set(gene for gene in GeneList))]

        if LineValues[6] == "CRISPR":
            flag += 1
            CRISPRcoordinates.append(count)
            if flag == 1:
                IslandID = str(LineValues[4] + '_' + '_'.join(LineValues[1].split('..')))
            #print(IslandID)
            for gene in GeneList[max(count - 4, 0):count-1]:
                UniqueSet.add(gene)
            IslandIDwithGenes[IslandID] = [list(UniqueSet), len(set(gene for gene in GeneList))]
            if flag > 4:
                print(flag, IslandID)
        if LineValues[6] == "Unknown":
            GeneList.append(LineValues[11])
            if (flag == 1 and len(UniqueSet) < min(6, 2 + CRISPRcoordinates[0])) \
                    or (flag == 2 and len(UniqueSet) < min(12, min(6, 2 + CRISPRcoordinates[0]) +
                        min(6, CRISPRcoordinates[1] - CRISPRcoordinates[0] - 1))) \
                    or ((flag == 3 and len(UniqueSet) < min(18, min(6, 2 + CRISPRcoordinates[0]) +
                        min(6, CRISPRcoordinates[1] - CRISPRcoordinates[0] - 1) +
                        min(6, CRISPRcoordinates[2] - CRISPRcoordinates[1] - 1)))) \
                    or (flag == 4 and len(UniqueSet) < min(24, min(6, 2 + CRISPRcoordinates[0]) +
                        min(6, CRISPRcoordinates[1] - CRISPRcoordinates[0] - 1) +
                        min(6, CRISPRcoordinates[2] - CRISPRcoordinates[1] - 1) +
                        min(6, CRISPRcoordinates[3] - CRISPRcoordinates[2] - 1))) \
                    or (flag == 5 and len(UniqueSet) < min(36, min(6, 2 + CRISPRcoordinates[0]) +
                        min(6, CRISPRcoordinates[1] - CRISPRcoordinates[0] - 1) +
                        min(6, CRISPRcoordinates[2] - CRISPRcoordinates[1] - 1) +
                        min(6, CRISPRcoordinates[3] - CRISPRcoordinates[2] - 1) +
                        min(6, CRISPRcoordinates[4] - CRISPRcoordinates[3] - 1))) :

                UniqueSet.add(LineValues[11])
            if IslandID in IslandIDwithGenes:
                IslandIDwithGenes[IslandID] = [list(UniqueSet), len(set(gene for gene in GeneList))]
    # print('UniqueSet:   ', UniqueSet, '     genelist: ', GeneList)
    # print(IslandIDwithGenes)

with open(IslandIDwithCRISPRNeighborsFileName, "w") as IslandsInfoFile:
    for ID in IslandIDwithGenes:
        IslandsInfoFile.write(ID + '\t' + Unlist(IslandIDwithGenes[ID][0]) + '\t'
                              + str(IslandIDwithGenes[ID][1]) + '\n')



CRISPRinfoByID = {}
for line in open(IdentifiedArraysInfoFileName,"r"):
    LineValues = line[:-1].split("\t")
    CRISPRinfoByID[str(LineValues[1].replace('"', '') + '_' + LineValues[2]+ '_' + LineValues[3])] = str('\t' + LineValues[4] + '\t' + LineValues[5].replace('"', '') + '\t' +
                                                                                        LineValues[10].replace('"', '') + '\t')
count1 = 0
ReplacingID = {}
with open( IslandIDwithCRISPRNeighborsFixedFileName, "w") as File:
    for line in open(IslandIDwithCRISPRNeighborsFileName, "r"):
        count1 += 1
        LineValues = line[:-1].split('\t')
        GeneList = ','.join(LineValues[1][1:].split(' '))
        ID = LineValues[0]
        Genome = LineValues[0].split('_')[0]
        for key in CRISPRinfoByID:
            if re.match(Genome, key.split('_')[0]) and ID.split('_')[1] == key.split('_')[1] and ID.split('_')[2] == key.split('_')[2]:  # Because in the Islands.ann_clust IDs are like "CP00001" while in initial table they are "CP00001.1"
                ReplacingID[ID] = key
                ID = key

        if ID in CRISPRinfoByID:
            CRISPRinfoByID[ID] += str(GeneList)
            File.write(ID + '\t' + str(CRISPRinfoByID[ID]) + '\n')
        print(count1)
print(len(CRISPRinfoByID))