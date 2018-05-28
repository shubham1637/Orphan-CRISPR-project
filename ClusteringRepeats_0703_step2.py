import math
import datetime
import numpy as np
import itertools
from scipy.cluster.hierarchy import dendrogram, linkage, average, fcluster, fclusterdata, to_tree


NonRedundantSequencesFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Nonredundant_repeats.txt"
RepeatsAssignedFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/RepeatsIdentical.txt"
RepeatsSimilaritiesFileName = "/panfs/pan1/orphancrispr/AllRepeatsVsAllRepeats_fix.hits"
AllArraysInfoFileName ="/home/utkinai2/Project1/merged_all.csv"
FilteredRepeatsSimilaritiesFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/AllRepeatsVsAllRepeats_fix_filtered_Dist.hits"


RepeatLengthsByID = {}
for line in open(AllArraysInfoFileName,"r"):
    LineValues = line[:-1].split(",")
    RepeatLengthsByID[str(LineValues[2] + '_' + LineValues[3] + '_' + LineValues[4])] = len(LineValues[6])

NonRedundantSet = set()
for line in open(NonRedundantSequencesFileName, "r"):
    NonRedundantSet.add(line[:-1])

UniqueSet = set()
print(datetime.datetime.now())
with open (FilteredRepeatsSimilaritiesFileName, "w") as FilteredFile:
    for line in open(RepeatsSimilaritiesFileName, "r"):
        LineValues = line[:-1].split('\t')
        Name1 = LineValues[0]
        Name2 = LineValues[1]
        if (Name1, Name2) in UniqueSet:
            continue

        if Name1 in NonRedundantSet and Name2 in NonRedundantSet:
            Coverage = int(LineValues[2]) / int(max(RepeatLengthsByID[Name2], RepeatLengthsByID[Name1]))
            if Coverage > 0.8:
                Score = int(LineValues[3]) * min(RepeatLengthsByID[Name1], RepeatLengthsByID[Name2]) / (
                    max(RepeatLengthsByID[Name2], RepeatLengthsByID[Name1]) * 2 * int(LineValues[2]))
                if Score != 1:
                    Dist = - math.log(Score)
                    FilteredFile.write(Name1 + '\t' + Name2 + '\t' + str(Dist) + '\n')
                    UniqueSet.add((Name1, Name2))
                else:
                    if (Name1,Name2) not in UniqueSet:
                        Dist = 0.0
                        FilteredFile.write(Name1 + '\t' + Name2 + '\t' + str(Dist) + '\n')
                        UniqueSet.add((Name1, Name2))
                    else:
                        continue
                    # if (Name1, Name2) in DictWithDoubledScores:
                    #     DictWithDoubledScores[Name1, Name2] = max(int(DictWithDoubledScores[Name1, Name2]), Dist)
                    # else:
                    #     DictWithDoubledScores[Name1, Name2] = max(0, Dist)

            # else:
            #     if (Name1, Name2) not in UniqueSet:
            #         Dist = 1.0
            #         FilteredFile.write(Name1 + '\t' + Name2 + '\t' + str(Dist) + '\n')
            #         UniqueSet.add((Name1, Name2))
            #     else:
            #         continue
print(datetime.datetime.now())






