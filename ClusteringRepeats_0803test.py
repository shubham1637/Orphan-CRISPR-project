import datetime
import numpy as np
import itertools
import scipy
from scipy.cluster.hierarchy import dendrogram, linkage, average, fcluster, fclusterdata, to_tree
import math
import io

NonRedundantSequencesFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Nonredundant_repeats_test.txt"
RepeatsAssignedFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/RepeatsIdentical_test.txt"
RepeatsSimilaritiesFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/TEST.hits"
AllArraysInfoFileName ="/home/utkinai2/Project1/merged_all.csv"
FilteredRepeatsSimilaritiesFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/AllRepeatsVsAllRepeats_test_filtered.hits"
NonRedundantRepeatsWithNumbersFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Nonredundant_repeats_with_numbers_test.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName100 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Test_NonredundantRepeats_with_ClusterLabels_1.0.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName90 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Test_NonredundantRepeats_with_ClusterLabels_0.9.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName80 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Test_NonredundantRepeats_with_ClusterLabels_0.8.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName70 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Test_NonredundantRepeats_with_ClusterLabels_0.7.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName50 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Test_NonredundantRepeats_with_ClusterLabels_0.5.txt"

# IDwithLabelsFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/All_SS_withClusterLabels_forRepresentatives_0.9.txt"
# AllIDwithLabelsFileName = "/home/utkinai2/Project1/ClusteringRepeats/All_SS_withClusterLabels_all_0.9.txt"

RepeatLengthsByID = {}
for line in open(AllArraysInfoFileName,"r"):
    LineValues = line[:-1].split(",")
    RepeatLengthsByID[str(LineValues[2] + '_' + LineValues[3] + '_' + LineValues[4])] = len(LineValues[6])

def ReverseComplement(seq):
   complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
   reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
   return reverse_complement

print(datetime.datetime.now())

# count = 0
# SequenceToClustDict = dict()
# ID = ""
# for Line in open("/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/AllRepeats.fna"):
#     count += 1
#     if count > 500:
#         break
#     if Line[0] == ">":
#        ID = Line[1:-1]
#     else:
#        Repeat = Line[:-1]
#        ReverseComplementRepeat = ReverseComplement(Repeat)
#
#        if Repeat in SequenceToClustDict:
#            SequenceToClustDict[Repeat].append(ID)
#        elif ReverseComplementRepeat in SequenceToClustDict:
#            SequenceToClustDict[ReverseComplementRepeat].append(ID)
#        else:
#            SequenceToClustDict[Repeat] = [ID]
#
# #creating
# with open(RepeatsAssignedFileName, "w") as NonRedundantFile:
#     for Repeat in SequenceToClustDict:
#         NonRedundantFile.write(Repeat + '\t' + ' '.join(SequenceToClustDict[Repeat]) + '\n')

# NonRedundantSet = set()
# for line in open(RepeatsAssignedFileName,"r"):
#     RepresentativeID = line[:-1].split('\t')[1].split(' ')[0]
#     NonRedundantSet.add(RepresentativeID)
#
# with open(NonRedundantSequencesFileName ,"w") as f:
#     for case in NonRedundantSet:
#         f.write(case + '\n')

# здесь вручную взяты 6 репитов, ниже их сиквенсы

NonRedundantSet = ['CP003287.1_14346_14672','CP003287.1_39917_40243','CP003659.1_57665_59539','CP003659.1_3250545_3251389','CP003659.1_4752371_4758593' ,'CP003659.1_5001036_5003639']
# >CP003287.1_14346_14672
# ATTGCAATTTCACTTACTCCCTATTAGGGATTGAAAC
# >CP003287.1_39917_40243
# GTTTCAATCCCTAATAGGGAGTAAGTGAAATTGCAAT
# >CP003659.1_57665_59539
# ATTGCAATTTTCATTAATCCCTATTAGGGATTGAAAC
# >CP003659.1_3250545_3251389
# CTTTCAACCCGCCCCACTCCGGGAGGGGTGTTGAAAC
# >CP003659.1_4752371_4758593
# GTTTCAGTCCCCTTGCGGGGATTATTTATTTGGAAAC
# >CP003659.1_5001036_5003639
# GTTTCAATCCCTAATAGGGATTATTTGAAATTTCAAC


# костыль в виде UniqueSet , чтобы избежать повторяющихся в .hits файле пар
DictWithDoubledScores = {}
UniqueSet = set()
print(datetime.datetime.now())
with open (FilteredRepeatsSimilaritiesFileName, "w") as FilteredFile:
    for line in open(RepeatsSimilaritiesFileName, "r"):
        LineValues = line[:-1].split('\t')
        Name1 = LineValues[0]
        Name2 = LineValues[1]

        if Name1 in NonRedundantSet and Name2 in NonRedundantSet:
            Coverage = int(LineValues[2]) / int(max(RepeatLengthsByID[Name2], RepeatLengthsByID[Name1]))
            if Coverage > 0.8:
                Score = int(LineValues[3]) * min(RepeatLengthsByID[Name1], RepeatLengthsByID[Name2]) / (
                    max(RepeatLengthsByID[Name2], RepeatLengthsByID[Name1]) * 2 * int(LineValues[2]))
                if Score != 1:
                    Dist = - math.log2(Score)
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
                #     DictWithDoubledScores[Name1, Name2] = min(int(DictWithDoubledScores[Name1, Name2]),Dist)
                # else:
                #     DictWithDoubledScores[Name1, Name2] =  Dist

print(datetime.datetime.now())

print(NonRedundantSet)
print(len(NonRedundantSet))
count = 0


# это делается, чтобы заполнять далее нужную ячейку эррэя: например, ID1 --> 1 и IDblabla --> 15,
# т.е. ячейка пары ID1:IDblabla будет array[1, 15], в которой будет лежать distance между парой
NonRedundantSetWithNumbers = {}
with open(NonRedundantRepeatsWithNumbersFileName, "w") as NonRedundantRepeatsWithNumbersFile:
    for Name in NonRedundantSet:
        NonRedundantSetWithNumbers[Name] = count
        NonRedundantRepeatsWithNumbersFile.write(Name + '\t' + str(count) +'\n')
        count += 1
print(NonRedundantSetWithNumbers)
print(datetime.datetime.now())

count1=0
DistanceMatrix = np.full((len(NonRedundantSet), len(NonRedundantSet)), 1.0)



for line in open(FilteredRepeatsSimilaritiesFileName, "r"):
    count1 += 1
    LineValues = line[:-1].split('\t')
    Name1 = LineValues[0]
    Name2 = LineValues[1]
    DistanceMatrix[NonRedundantSetWithNumbers[Name1], NonRedundantSetWithNumbers[Name2]] = LineValues[2]
print(DistanceMatrix)

np.savetxt('/panfs/pan1/orphancrispr/ClusteringRepeats/0703/DistanceMatrixTest.txt', DistanceMatrix, delimiter = ',')

#DistanceMatrix = np.reshape(DistanceMatrix, (36,1))   - это попытка сделать condensed matrix , т.е. просто вектор
#DistanceMatrix = np.ndarray.flatten(DistanceMatrix)    - ещё одна

# scipy.spatial.distance.squareform(DistanceMatrix)   - эта функция обещает сделать всё, что нужно, но это обман.

#print(DistanceArray)
print(DistanceMatrix)

print(datetime.datetime.now(), ' Distance Matrix is created')

#DistanceMatrixUpperDiag = np.triu(DistanceMatrix, k=0)
#print(DistanceMatrixUpperDiag)

print(datetime.datetime.now(), ' Start of hierarchical clustering with distance = 0')
SimilarityMatrixArray = linkage(DistanceMatrix, method = 'single')
FclusterArray = fcluster(SimilarityMatrixArray, 0.0, criterion='distance', depth=2, R=None, monocrit=None)
print(datetime.datetime.now(), ' End of hierarchical clustering with distance = 0')

with open(NonRedundantRepeatsIDwithClusterLabelsFileName100, "w") as File:
    i=0
    for name in NonRedundantSet:
       File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
       i += 1

print(datetime.datetime.now(), ' Start of hierarchical clustering with distance = 0.1')
SimilarityMatrixArray = linkage(DistanceMatrix, method = 'single')
FclusterArray = fcluster(SimilarityMatrixArray, 0.1, criterion='distance', depth=2, R=None, monocrit=None)
print(datetime.datetime.now(), ' End of hierarchical clustering with distance = 0.1')

with open(NonRedundantRepeatsIDwithClusterLabelsFileName90, "w") as File:
    i=0
    for name in NonRedundantSet:
       File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
       i += 1

print(datetime.datetime.now(), ' Start of hierarchical clustering with distance = 0.2')
SimilarityMatrixArray = linkage(DistanceMatrix, method = 'single')
FclusterArray = fcluster(SimilarityMatrixArray, 0.2, criterion='distance', depth=2, R=None, monocrit=None)
print(datetime.datetime.now(), ' End of hierarchical clustering with distance = 0.2')

with open(NonRedundantRepeatsIDwithClusterLabelsFileName80, "w") as File:
    i=0
    for name in NonRedundantSet:
       File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
       i += 1

print(datetime.datetime.now(), ' Start of hierarchical clustering with distance = 0.3')
SimilarityMatrixArray = linkage(DistanceMatrix, method = 'single')
FclusterArray = fcluster(SimilarityMatrixArray, 0.5, criterion='distance', depth=2, R=None, monocrit=None)
print(datetime.datetime.now(), ' End of hierarchical clustering with distance = 0.3')

with open(NonRedundantRepeatsIDwithClusterLabelsFileName70, "w") as File:
    i=0
    for name in NonRedundantSet:
       File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
       i += 1

print(datetime.datetime.now(), ' Start of hierarchical clustering with distance = 0.5')
SimilarityMatrixArray = linkage(DistanceMatrix, method = 'single')
FclusterArray = fcluster(SimilarityMatrixArray, 1, criterion='distance', depth=2, R=None, monocrit=None)
print(datetime.datetime.now(), ' End of hierarchical clustering with distance = 0.5')

with open(NonRedundantRepeatsIDwithClusterLabelsFileName50, "w") as File:
    i=0
    for name in NonRedundantSet:
       File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
       i += 1

CRISPRinfoByID = {}
for line in open(AllArraysInfoFileName,"r"):
    LineValues = line[:-1].split(",")
    CRISPRinfoByID[str(LineValues[2] + '_' + LineValues[3]+ '_' + LineValues[4])] = str(LineValues[5] + '\t' + LineValues[6] + '\t' +
                                                                                        LineValues[10] + '\t' + LineValues[11] + '\t' +
                                                                                        LineValues[15])

ClustersDict = {}
Similarity = [0.5,0.7,0.8,0.9,1.0]
for threshold in Similarity:
    with io.open('/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Clusters_ofNonRedundantRepeats_test' + str(threshold) + '.txt', "w") as File:
        print(threshold)
        for line in io.open("/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Test_NonredundantRepeats_with_ClusterLabels_" + str(threshold) + ".txt" , "r"):
            LineValues  = line[:-1].split('\t')
            if LineValues[0] in CRISPRinfoByID:
                File.write(str(LineValues[1]) + '\t' + LineValues[0] + '\t' + CRISPRinfoByID[LineValues[0]] + '\n')