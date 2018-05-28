import math
import datetime
import numpy as np
import pandas as pd
import itertools

SpacerIDWithNumberFileName = "/panfs/pan1/orphancrispr/SpacersFromArraysCluster/OrphanSpacerID_withNumber.txt"
SpacersSimilaritiesFileName = "/panfs/pan1/orphancrispr/SpacersFromArraysCluster/OrphanSpacerHits_filtered90-90.hits"
# SpacersSimilaritiesWGFileName = "/panfs/pan1/orphancrispr/SpacersFromArraysCluster/IdentifiedSpacerHits_filtered90-90_wholegenomes.hits"
SpacersWithDistanceFileName = "/panfs/pan1/orphancrispr/SpacersFromArraysCluster/OrphanSpacerHits_withDistance_2003.txt"
# SpacersWithDistanceFileName = "/panfs/pan1/orphancrispr/SpacersFromArraysCluster/TestFile1.txt"
AssignedToNonRedundantFileName = "/panfs/pan1/orphancrispr/SpacersFromArraysCluster/NonRedundantSpacers_withAssigned_Orphans_2803.txt"
# WholeGenomesFileName = "/panfs/pan1/orphancrispr/Prok1603.pp.txt"

# run this part only once to create a file
print(datetime.datetime.now())
# with open(SpacersSimilaritiesFileName, "w") as FilteredSpacerHitsFile:
#     for line in open("/panfs/pan1/orphancrispr/SpacersFromArraysCluster/Spacers_from_Orphans_2003.hits" ,"r"):
#         if line[0] == "#":
#             continue
#         LineValues = line[:-1].split('\t')
#         if float(LineValues[6]) >= 90.0 and float(LineValues[7]) >= 90.0:
#             #print('YES' , ID[0], ' ', identity, '% ', coverage, '%')
#             FilteredSpacerHitsFile.write(line)
#
# # run this part only once to create a file
# print(datetime.datetime.now(), ' SpacerHits filtered')
# WholeGenomes = set()
# for line in open(WholeGenomesFileName, "r"):
#     WholeGenomes.add(line[:-1].split('\t')[1])
#
# with open(SpacersSimilaritiesWGFileName, "w") as File:
#     for line in open(SpacersSimilaritiesFileName, "r"):
#         if line[:-1].split('\t')[0].split('_')[0] in WholeGenomes:
#             File.write(line)
# print(datetime.datetime.now(), ' SpacerHits filtered by wg')
#
# with open(SpacersWithDistanceFileName, "w") as File:
#     for line in open(SpacersSimilaritiesWGFileName, "r"):
#         LineValues = line[:-1].split('\t')
#         Name1 = LineValues[0]
#         Name2 = LineValues[1]
#         Score = (float(LineValues[6])/100) * (float(LineValues[7])/100)
#         if Score != 1:
#             Dist = - math.log(Score)
#             File.write(Name1 + '\t' + Name2 + '\t' + str(Dist) + '\n')
#         else:
#             Dist = 0.0
#             File.write(Name1 + '\t' + Name2 + '\t' + str(Dist) + '\n')
#
#
# print(datetime.datetime.now(), 'Distances calculated')

AllSpacersID = set()
NonRedundantSetSpacers = {}
NonRedundantSet = set()
AssignedSet = set()
Set = set()
count = 0
SelfMatch = set()
for line in open(SpacersWithDistanceFileName, "r"):
    LineValues = line[:-1].split('\t')
    AllSpacersID.add(LineValues[0])

    if float(LineValues[2]) == 0.0:
            if LineValues[0] == LineValues[1]:
                if LineValues[0] not in AssignedSet:
                    SelfMatch.add(LineValues[0])
                    continue
            if LineValues[0] in AssignedSet and LineValues[0] not in NonRedundantSet:
                if LineValues[1] not in NonRedundantSet and LineValues[1] not in AssignedSet:
                    for key in NonRedundantSetSpacers:
                        if LineValues[0] in NonRedundantSetSpacers[key]:
                            NonRedundantSetSpacers[key] += [LineValues[1]]
                            AssignedSet.add(LineValues[1])

            if LineValues[0] in NonRedundantSet and LineValues[0] not in AssignedSet and LineValues[1] not in AssignedSet \
                and LineValues[1] not in NonRedundantSet:
                if LineValues[0] in NonRedundantSetSpacers:
                    NonRedundantSetSpacers[LineValues[0]] += [LineValues[1]]
                else:
                    NonRedundantSetSpacers[LineValues[0]] = [LineValues[1]]
                AssignedSet.add(LineValues[1])

            if LineValues[0] not in NonRedundantSetSpacers and LineValues[0] not in AssignedSet:

                if LineValues[1] in AssignedSet and LineValues[1] not in NonRedundantSet:
                    for key in NonRedundantSetSpacers:
                        if LineValues[1] in NonRedundantSetSpacers[key]:
                            NonRedundantSetSpacers[key] += [LineValues[0]]
                    AssignedSet.add(LineValues[0])

                if LineValues[1] in NonRedundantSet and LineValues[1] not in AssignedSet:
                    if LineValues[1] in NonRedundantSetSpacers:
                        NonRedundantSetSpacers[LineValues[1]] += [LineValues[0]]
                    else:
                        NonRedundantSetSpacers[LineValues[1]] = [LineValues[0]]
                    AssignedSet.add(LineValues[0])

                if LineValues[1] not in NonRedundantSet and LineValues[1] not in AssignedSet:
                    NonRedundantSet.add(LineValues[0])
                    NonRedundantSetSpacers[LineValues[0]] = [LineValues[1]]
                    AssignedSet.add(LineValues[1])

for line in open(SpacersWithDistanceFileName, "r"):
    LineValues = line[:-1].split('\t')
    if float(LineValues[2]) != 0.0:
        if LineValues[0] not in AssignedSet:
            NonRedundantSet.add(LineValues[0])
        if LineValues[1] not in AssignedSet:
            NonRedundantSet.add(LineValues[1])

for spacer in SelfMatch:
    if spacer not in AssignedSet:
        NonRedundantSet.add(spacer)

SelfMatchOnly = set(SelfMatch - NonRedundantSet)

# NonRedundantSetSpacers1 = {}
# for key in NonRedundantSetSpacers:
#     for value in NonRedundantSetSpacers[key]:
#         if value in NonRedundantSet:
#             continue
#         else:
#             if key not in NonRedundantSetSpacers1:
#                 NonRedundantSetSpacers1[key] = [value]
#             else:
#                 NonRedundantSetSpacers1[key] += [value]
# Length = 0
# for key in NonRedundantSetSpacers1:
#     Length += len(NonRedundantSetSpacers1[key])
print('Self match only ', len(SelfMatchOnly))
print('Dict length = ', len(NonRedundantSetSpacers))
print(' NonRedundant Set ' ,len(NonRedundantSet))
print(' Assigned Set ', len(AssignedSet))
print('intersection between NonRed and Assigned ', len(list(AssignedSet & NonRedundantSet)))

SpacerFasta = {}
with open('/panfs/pan1/orphancrispr/SpacersFromArraysCluster/Spacers_from_Orphans.txt', "r") as f:
    for line1, line2 in itertools.zip_longest(* [f] * 2):
         SpacerFasta[line1[:-1]] = line2[:-1]

with open('/panfs/pan1/orphancrispr/SpacersFromArraysCluster/OrphanSpacers_NonRedundantSetOnly_2803.txt', "w") as File:
    for nonred in NonRedundantSet:
        File.write(nonred + '\n' + SpacerFasta[nonred] + '\n')



with open(AssignedToNonRedundantFileName, "w") as File:
    for spacer in NonRedundantSetSpacers:
        File.write(spacer + '\t' + str((NonRedundantSetSpacers[spacer])) +'\n')
    for spacer in NonRedundantSet:
        if spacer not in NonRedundantSetSpacers:
            File.write(spacer + '\n')

# count = 0
# CountedSpacersForNumpy = {}
# NumberIntoID = {}
# with open(SpacerIDWithNumberFileName, "w") as File:
#     for ID in NonRedundantSetSpacers:
#         CountedSpacersForNumpy[ID] = count
#         NumberIntoID[count] = ID
#         File.write(ID + '\t' + str(count) +'\n')
#         count += 1
#
# print(count)
print(datetime.datetime.now(), '  all spacers: ',  len(AllSpacersID))
#
# DistanceMatrix = np.ones((len(NonRedundantSet), len(NonRedundantSet)))
# #np.full((2, 2), 10)
#
# for line in open(SpacersWithDistanceFileName , "r"):
#     LineValues = line[:-1].split('\t')
#     ID1 = LineValues[0]
#     ID2 = LineValues[1]
#     if ID1 in CountedSpacersForNumpy and ID2 in CountedSpacersForNumpy:
#         DistanceMatrix[CountedSpacersForNumpy[ID1], CountedSpacersForNumpy[ID2]] = LineValues[2]
#
# print(datetime.datetime.now(), ' Distance Matrix is created')
#
# df = pd.DataFrame(DistanceMatrix)
# df = df.rename(index=NumberIntoID, columns=NumberIntoID)
# df.to_csv("/panfs/pan1/orphancrispr/SpacersFromArraysCluster/DistanceMatrix_spacers_Identified_.csv")
#
# print(datetime.datetime.now())