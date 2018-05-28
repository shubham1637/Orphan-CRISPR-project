import re

OrphanArraysInfoFileName ="/home/utkinai2/Project1/Orphans_finalfiltered.txt"
IslandIDwithCRISPRNeighborsFileName = "/panfs/pan1/orphancrispr/IslandsCluster/Islands_pfam_3NeighborsSet.txt"
IslandIDwithCRISPRNeighborsFixedFileName = "/panfs/pan1/orphancrispr/IslandsCluster/Islands_pfam_3NeighborsSet_fixedID.txt"

CRISPRinfoByID = {}
for line in open(OrphanArraysInfoFileName,"r"):
    LineValues = line[:-1].split(" ")
    CRISPRinfoByID[str(LineValues[2].replace('"', '') + '_' + LineValues[3]+ '_' + LineValues[4])] = str('\t' + LineValues[5] + '\t' + LineValues[6].replace('"', '') + '\t' +
                                                                                        LineValues[11].replace('"', '') + '\t')
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

# for line in open(IslandIDwithCRISPRNeighborsFixedFileName, "r"):
#     ID = line[:-1].split('\t')[0]
#     GeneList = line[:-1].split('\t')[4]
#     if ID in CRISPRinfoByID:
#         CRISPRinfoByID[ID] += str(GeneList)
#
#
# Threshold = [0, 0.2, 0.5, 0.75,  1, 1.25, 1.5, 1.6, 1.65, 1.7, 1.75 , 1.8, 1.9, 2, 2.25, 2.5, 3, 3.5, 4, 4.5, 5]
# for threshold in Threshold:
#     with open("/panfs/pan1/orphancrispr/IslandsCluster/Output_pfamclustered_h_" + str(threshold) + ".txt","w") as File:
#         count = 0
#         for line in open("/panfs/pan1/orphancrispr/IslandsCluster/Pfamgroups_h_" + str(threshold) + ".txt","r"):
#             print(threshold)
#             count +=1
#             if count == 1:
#                 continue
#
#             LineValues = line[:-1].split(' ')
#             ID = LineValues[0].strip('"')
#             cluster = LineValues[1]
#             File.write(str(cluster) + '\t' + '\t'.join(ID.split('_')) + CRISPRinfoByID[ID] + '\n')
#
