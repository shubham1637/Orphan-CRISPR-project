import math
import numpy as np
import pandas as pd
import datetime

# IslandsAnnotationFileName = "/panfs/pan1/orphancrispr/IslandsCluster/Islands_pfam.ann_clust"
# IslandIDwithCRISPRNeighborsFileName = "/panfs/pan1/orphancrispr/IslandsCluster/Islands_pfam_3NeighborsSet_fixedID.txt"
# PairwiseIslandDistancesFileName = "/panfs/pan1/orphancrispr/IslandsCluster/PairwiseIslandDistances_pfam.txt"
# IslandIDWithNumberFileName = "/panfs/pan1/orphancrispr/IslandsCluster/IslandsID_pfam_withNumber.txt"
IslandsAnnotationFileName = "/panfs/pan1/prokdata/CRISPRicity/JointClusters/MergedIslands.tsv"
IslandIDwithCRISPRNeighborsFileName = "/panfs/pan1/orphancrispr/IslandsCluster/Islands_Identified_3NeighborsSet_fixedID.txt"
PairwiseIslandDistancesFileName = "/panfs/pan1/orphancrispr/IslandsCluster/PairwiseIslandDistances_Identified.txt"
IslandIDWithNumberFileName = "/panfs/pan1/orphancrispr/IslandsCluster/IslandsID_Identified_withNumber.txt"

def intersect(List1, List2):
    return list(set(List1) & set(List2))


AllIDs = set()
CountedIDsForNumpy = {}
NumberIntoID = {}
count = 0
with open(IslandIDWithNumberFileName, "w") as IslandIDWithNumberFile:
    for line in open(IslandIDwithCRISPRNeighborsFileName, "r"):
        ID = line[:-1].split('\t')[0]
        AllIDs.add(ID)
        CountedIDsForNumpy[ID] = count
        NumberIntoID[count] = ID
        IslandIDWithNumberFile.write(ID + '\t' + str(count) +'\n')
        count += 1

IDandGenesDict = {}
for line in open(IslandIDwithCRISPRNeighborsFileName, "r"):
    LineValues = line[:-1].split('\t')
    IDandGenesDict[LineValues[0]] = [LineValues[4].split(',')]


print(datetime.datetime.now())
with open(PairwiseIslandDistancesFileName, "w") as DistanceFile:
    for ID1 in AllIDs:
        for ID2 in AllIDs:
            #print(ID1, ' ', ID2, ' ', intersect(IDandGenesDict[ID1][0], IDandGenesDict[ID2][0]))
            Score = len(intersect(IDandGenesDict[ID1][0], IDandGenesDict[ID2][0])) / \
                    len(list(set().union(IDandGenesDict[ID1][0],IDandGenesDict[ID2][0])))

            if Score == 1.0:
                Dist = 0.0
            if Score == 0:
                Dist = 5.0
            else:
                if Score != 1.0:
                    Dist = - math.log(Score)
            #print(ID1, '  ', ID2, ' ', Dist)
            DistanceFile.write(ID1 + '\t' + ID2 + '\t' + str(Dist) + '\n')
print(datetime.datetime.now(), ' Pairwise distances calculated')

DistanceMatrix = np.full((len(AllIDs), len(AllIDs)), 5.0)

for line in open(PairwiseIslandDistancesFileName, "r"):
    LineValues = line[:-1].split('\t')
    ID1 = LineValues[0]
    ID2 = LineValues[1]
    DistanceMatrix[CountedIDsForNumpy[ID1], CountedIDsForNumpy[ID2]] = LineValues[2]

print(datetime.datetime.now(), ' Distance Matrix is created')


df = pd.DataFrame(DistanceMatrix)
df = df.rename(index=NumberIntoID, columns=NumberIntoID)
df.to_csv("/panfs/pan1/orphancrispr/IslandsCluster/DistanceMatrix_IdentifiedIslands.csv")
