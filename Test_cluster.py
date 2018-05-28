import datetime
import pandas as pd
import numpy as np
import itertools
from scipy.cluster.hierarchy import dendrogram, linkage, average, fcluster, fclusterdata, to_tree
from scipy.cluster.hierarchy import cophenet


def ReverseComplement(seq):
   complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
   reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
   return reverse_complement

RepeatsSimilarities3FileName = "//frosty/utkinai2/Project1/AllRepeatsVsAllRepeats.hits"
#RepeatsSimilaritiesFileName = "/panfs/pan1/orphancrispr/AllRepeatsVsAllRepeats_fix.hits"
RepeatsSimilaritiesFileName = "//frosty/utkinai2/Project1/ClusteringRepeats/ClusterableSet_UPGMA/ClusterableSetAllvsAll.hits"
RepeatsSimilarities10000FileName = "/panfs/pan1/orphancrispr/Test/17000000Test.hits"
NonRedundatRepeatIDFile = "//frosty/utkinai2/Project1/ClusteringRepeats/ClusterableSet_UPGMA/NonRedundantRepeatIDClusterableSet.txt" # "/panfs/pan1/orphancrispr/Test/NonRedundantSet_17m.txt"
PairwiseSimilaritiesDictFileName = "//frosty/utkinai2/Project1/ClusteringRepeats/ClusterableSet_UPGMA/PairwiseSim_ClusterableSet.txt" # "/panfs/pan1/orphancrispr/Test/Pairwise_similaritiesDict_17m.txt"
AssignedDictFileName =    "//frosty/utkinai2/Project1/ClusteringRepeats/ClusterableSet_UPGMA/AssignedIDs_ClusterableSet.txt"   # "/panfs/pan1/orphancrispr/Test/AssignedIDs_17m.txt"
AllArraysInfoFileName ="//frosty/utkinai2/Project1/merged_all.csv"
ClusterableSetInfoFileName = "//frosty/utkinai2/Project1/ClusteringRepeats/ClusterableSet_UPGMA/ClusterableSet_fasta.txt"
RepresentativeIDwithLabelsFileName = "//frosty/utkinai2/Project1/ClusteringRepeats/ClusterableSet_UPGMA/ClusterableSet_forRepresentative.txt"
AllIDwithLabelsFileName = "//frosty/utkinai2/Project1/ClusteringRepeats/ClusterableSet_UPGMA/ClusterableSet_WithClusterLabels_all_0.99_dist.txt"
NameListInDigitsFileName = "//frosty/utkinai2/Project1/ClusteringRepeats/ClusterableSet_UPGMA/NonRedundatRepeatsIDinDigits.txt"

RepeatLengthsByID_3 = {'A': 30, 'B': 27, 'C': 27, 'D':27}
RepeatLengthsByID = {}


with open(ClusterableSetInfoFileName, "r") as f:
    for line1, line2 in itertools.zip_longest(*[f] * 2):
        RepeatLengthsByID[line1[:-1].strip('>')] = len(line2[:-1])
# for line in open(AllArraysInfoFileName,"r"):
#     LineValues = line[:-1].split(",")
#     RepeatLengthsByID[str(LineValues[2] + '_' + LineValues[3] + '_' + LineValues[4])] = len(LineValues[6])

# NameList0 = []
# for line in open(AllFileName, "r"):
#     LineValues = line[:-1].split('\t')
#     NameList0.append(str(LineValues[2] +'_' + LineValues[3] +'_' + LineValues[4]))
# print(len(NameList0))

1234
NonRedundantSet = set()
AssignedDict = {}
AssignedDictValues = []
count = 0
PairwiseSimilaritiesDict = {}
# NonRedundantSet = set()
print(datetime.datetime.now())
for line in open(RepeatsSimilaritiesFileName, "r"):
    #try:
        count += 1
        # if count % 100000 == 0:
        #     print(count)
        LineValues = line[:-1].split('\t')
        Name1 = LineValues[0]
        Name2 = LineValues[1]
        if Name1 in AssignedDictValues:
            continue
        if Name1 == Name2:
            PairwiseSimilaritiesDict[Name1, Name2] = 1
            NonRedundantSet.add(Name1)
            continue
        Score = int(LineValues[3]) * min(RepeatLengthsByID[Name1], RepeatLengthsByID[Name2])/ (max(RepeatLengthsByID[Name2], RepeatLengthsByID[Name1]) * 2 * int(LineValues[2]))
        if Score == 1:
            if Name1 not in AssignedDict.keys():
                AssignedDict[Name1] = [Name2]
                AssignedDictValues.append(Name2)
            else:
                if [Name2] not in AssignedDict[Name1]:
                    AssignedDict[Name1] += [Name2]
                    AssignedDictValues.append(Name2)
        else:
            if [Name1] in AssignedDict.values():
                continue
            NonRedundantSet.add(Name1)
    # except:
    #     print(count)
for Name in AssignedDict.keys():
    NonRedundantSet.add(Name)

print(NonRedundantSet)

for line in open(RepeatsSimilaritiesFileName, "r"):
        count += 1
        LineValues = line[:-1].split('\t')
        Name1 = LineValues[0]
        Name2 = LineValues[1]
        if Name1 in NonRedundantSet and Name2 in NonRedundantSet:
            Score = int(LineValues[3]) * min(RepeatLengthsByID[Name1], RepeatLengthsByID[Name2])/ (max(RepeatLengthsByID[Name2], RepeatLengthsByID[Name1]) * 2 * int(LineValues[2]))
            PairwiseSimilaritiesDict[Name1, Name2] = Score




print('length of NonRedundantSet  ',len(NonRedundantSet))
print('length of Dictionary with cluster representatives ', len(AssignedDict))
print(datetime.datetime.now())
print('lengths of PairwiseSimilaritiesDict ', len(PairwiseSimilaritiesDict.keys()))
print(PairwiseSimilaritiesDict)

count1 = 0
NonRedundantSetInDigits = {}
with open(NameListInDigitsFileName, "w") as NameListInDigitsFile:
    for Name in NonRedundantSet:
        NonRedundantSetInDigits[Name] = count1
        NameListInDigitsFile.write(Name + '\t' + str(count1))
        count1 += 1
    print(count1)

# SimilarityDataFrame = pd.DataFrame(index=NonRedundantSet, columns=NonRedundantSet)
# SimilarityDataFrame = SimilarityDataFrame.fillna(0)

# Assigning scores from PairwiseSimilaritiesDict to the corresponding indices to the zero-matrix
ScoreMatrix = np.zeros((len(NonRedundantSet), len(NonRedundantSet)))
for (Name1, Name2) in PairwiseSimilaritiesDict:
    ScoreMatrix[NonRedundantSetInDigits[Name1], NonRedundantSetInDigits[Name2]] = PairwiseSimilaritiesDict[Name1, Name2]
    ScoreMatrix[NonRedundantSetInDigits[Name2], NonRedundantSetInDigits[Name2]] = PairwiseSimilaritiesDict[Name1, Name2]
#print(ScoreMatrix)
#print(SimilarityDataFrame)
print(datetime.datetime.now())

with open(PairwiseSimilaritiesDictFileName,"w" ) as f:
    for key in PairwiseSimilaritiesDict:
       f.write(str(key) + '\t' + str(PairwiseSimilaritiesDict[key]) + '\n')

with open(AssignedDictFileName, "w") as AssignedFile:
    for key in AssignedDict:
        AssignedFile.write(str(key) + '\t' + str(set(AssignedDict[key])) + '\n')


# UPGMA
SimilarityMatrixArray = linkage(ScoreMatrix, method = 'average')
# #Z = linkage(SimilarityDataFrame, method = 'average')
# print(SimilarityMatrixArray)
#print(dendrogram(SimilarityMatrixArray))
#print(to_tree(SimilarityMatrixArray))
FclusterArray = fcluster(SimilarityMatrixArray, 0.1, criterion='distance', depth=2, R=None, monocrit=None)
#FclusterArray = fclusterdata(ScoreMatrix, 0.2, criterion = 'distance')
print(FclusterArray)


with open(RepresentativeIDwithLabelsFileName, "w") as File:
    i=0
    for name in list(NonRedundantSet):
       File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
       i += 1

ClustersDict = {}
with open(AllIDwithLabelsFileName, "w") as File:
    for line in open(RepresentativeIDwithLabelsFileName, "r"):
        TempList = set()
        LineValues = line[:-1].split('\t')
        if LineValues[0] in AssignedDict:                            # if this ID is a representative for previosly formed clusters (score =1)
            for i in range(len(AssignedDict[LineValues[0]])):
                TempList.add(str(AssignedDict[LineValues[0]][i]))    # adding to the temporary list all IDs assigned to the representative one (all with Score = 1 with representative)
            print(TempList)
            TempList = list(TempList)
        if len(TempList) != 0:
            if LineValues[1] in ClustersDict:
                ClustersDict[LineValues[1]] += [LineValues[0],str([str(TempList[i]) for i in range(len(TempList))])[2:-2]]
            else:
                ClustersDict[LineValues[1]] = [LineValues[0], str([str(TempList[i]) for i in range(len(TempList))])[2:-2]]
        else:
            if LineValues[1] in ClustersDict:
                ClustersDict[LineValues[1]] += [LineValues[0]]
            else:
                ClustersDict[LineValues[1]] = [LineValues[0]]

    for cluster in sorted(ClustersDict):
        File.write('Cluster ' + str(cluster) + '\n' + str(ClustersDict[cluster]).replace('"','') + '\n')
print('well, now enjoy this shitty clustering results')
print(ClustersDict)