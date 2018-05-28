import datetime
import numpy as np
import pandas as pd
import math

# IslandIDwithCRISPRNeighborsFixedFileName = "/panfs/pan1/orphancrispr/IslandsCluster/Islands_pfam_3NeighborsSet_fixedID.txt"
# IslandIDwithCRISPRNeighborsFixedFileName = "/panfs/pan1/orphancrispr/IslandsCluster/Test_neighborhood.txt"
#ArrayWithSpacersFileName =  "/panfs/pan1/orphancrispr/SpacersFromArraysCluster/OrphanArrays_with_Spacers.txt"
#ArrayWithSpacersFileName =  "/panfs/pan1/orphancrispr/Test_spacers.txt"
 #ArrayWithClusteredSpacersFileName =  "/panfs/pan1/orphancrispr/SpacersFromArraysCluster/OrphanArrays_with_ClusteredSpacers_2003.txt"
ArrayWithClusteredSpacersFileName =  "/panfs/pan1/orphancrispr/SpacersFromArraysCluster/IdentifiedArrays_with_ClusteredSpacers.txt"
IslandIDWithNumberFileName = "/panfs/pan1/orphancrispr/IslandsCluster/IslandsID_spacers_withNumber_Identified_wg.txt"
PairwiseIslandDistancesAfterSpacersClustFileName = "/panfs/pan1/orphancrispr/IslandsCluster/PairwiseIslandDistances_ClusteredSpacers_Identified_wg.txt"


def intersect(List1, List2):
    return list(set(List1) & set(List2))

IDandSpacers = {}
AllIDs = set()
CountedIDsForNumpy = {}
NumberIntoID = {}
count = 0

with open(IslandIDWithNumberFileName, "w") as IslandIDWithNumberFile:
    for line in open(ArrayWithClusteredSpacersFileName, "r"):
        LineValues = line[:-1].split('\t')
        ID = LineValues[0]
        IDandSpacers[ID] = LineValues[1].split(',')
        AllIDs.add(ID)
        CountedIDsForNumpy[ID] = count
        NumberIntoID[count] = ID
        IslandIDWithNumberFile.write(ID + '\t' + str(count) + '\n')
        count += 1
print(len(IDandSpacers.keys()))
print(len(AllIDs))

print(datetime.datetime.now())
with open(PairwiseIslandDistancesAfterSpacersClustFileName, "w") as DistanceFile:
    for ID1 in AllIDs:
        for ID2 in AllIDs:
            Score = len(intersect(IDandSpacers[ID1], IDandSpacers[ID2])) / \
                    len(list(set().union(IDandSpacers[ID1], IDandSpacers[ID2])))
            # print(ID1, ' ', ID2, ' ', intersect(IDandSpacers[ID1], IDandSpacers[ID2]), ' ', Score)
            if Score == 1.0:
                Dist = 0.0
            if Score == 0:
                Dist = 5.0
            else:
                if Score != 1.0:
                    Dist = - math.log(Score)
            # print(ID1, '  ', ID2, ' ', Dist)
            DistanceFile.write(ID1 + '\t' + ID2 + '\t' + str(Dist) + '\n')
print(datetime.datetime.now(), ' Pairwise distances calculated')

DistanceMatrix = np.full((len(AllIDs), len(AllIDs)), 5.0)

for line in open(PairwiseIslandDistancesAfterSpacersClustFileName, "r"):
    LineValues = line[:-1].split('\t')
    ID1 = LineValues[0]
    ID2 = LineValues[1]
    DistanceMatrix[CountedIDsForNumpy[ID1], CountedIDsForNumpy[ID2]] = LineValues[2]

print(datetime.datetime.now(), ' Distance Matrix is created')


df = pd.DataFrame(DistanceMatrix)
df = df.rename(index=NumberIntoID, columns=NumberIntoID)
df.to_csv("/panfs/pan1/orphancrispr/IslandsCluster/DistanceMatrix_withClusteredSpacers_IdentifiedWG.csv")
