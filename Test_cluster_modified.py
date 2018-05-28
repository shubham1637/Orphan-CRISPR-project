import datetime
import pandas as pd
import numpy as np
#from matplotlib import pyplot as plt
import itertools
import scipy
from scipy.cluster.hierarchy import dendrogram, linkage, average, fcluster, fclusterdata, to_tree
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist

RepeatsSimilarities3FileName = "//frosty/utkinai2/Project1/AllRepeatsVsAllRepeats.hits"
RepeatsSimilaritiesFileName = "/panfs/pan1/orphancrispr/AllRepeatsVsAllRepeats_fix.hits"
#RepeatsSimilaritiesFileName = "//frosty/utkinai2/Project1/ClusteringRepeats/ClusterableSetAllvsAll.hits"
RepeatsSimilarities10000FileName = "/panfs/pan1/orphancrispr/Test/17000000Test.hits"
NonRedundatRepeatIDFile =  "/panfs/pan1/orphancrispr/ClusteringRepeats/NonRedundantSetRepeats_0603.txt"
PairwiseSimilaritiesDictFileName =  "/panfs/pan1/orphancrispr/ClusteringRepeats/Pairwise_similaritiesDict_0603.txt"
AssignedDictFileName =   "/panfs/pan1/orphancrispr/Test/AssignedIDs.txt"
AllArraysInfoFileName ="/home/utkinai2/Project1/merged_all.csv"
ClusterableSetInfoFileName = "//frosty/utkinai2/Project1/ClusteringRepeats/ClusterableSet_fasta.txt"
RepresentativeIDwithLabelsFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/ClusterableSet_forRepresentative_0603.txt"
AllIDwithLabelsFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/WithClusterLabels_all_0.8_dist_0603.txt"
NameListInDigitsFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/NonRedundatRepeatsIDinDigits_0603.txt"

RepeatLengthsByID_3 = {'A': 30, 'B': 27, 'C': 27, 'D':27}
RepeatLengthsByID = {}

# with open(ClusterableSetInfoFileName, "r") as f:
#     for line1, line2 in itertools.zip_longest(*[f] * 2):
#         RepeatLengthsByID[line1[:-1].strip('>')] = len(line2[:-1])
for line in open(AllArraysInfoFileName,"r"):
    LineValues = line[:-1].split(",")
    RepeatLengthsByID[str(LineValues[2] + '_' + LineValues[3] + '_' + LineValues[4])] = len(LineValues[6])

# NameList0 = []
# for line in open(AllFileName, "r"):
#     LineValues = line[:-1].split('\t')
#     NameList0.append(str(LineValues[2] +'_' + LineValues[3] +'_' + LineValues[4]))
# print(len(NameList0))

# def plot_tree(dendro, pos=None ):
#     icoord = scipy.array( dendro['icoord'] )
#     dcoord = scipy.array( dendro['dcoord'] )
#     color_list = scipy.array( dendro['color_list'] )
#     xmin, xmax = icoord.min(), icoord.max()
#     ymin, ymax = dcoord.min(), dcoord.max()
#     if pos:
#         icoord = icoord[pos]
#         dcoord = dcoord[pos]
#         color_list = color_list[pos]
#     for xs, ys, color in zip(icoord, dcoord, color_list):
#         plt.plot(xs, ys,  color)
#     plt.xlim( xmin-10, xmax + 0.1*abs(xmax) )
#     plt.ylim( ymin, ymax + 0.1*abs(ymax) )
#     plt.show()

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
        if count % 100000 == 0:
            print(count)
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

NonRedundantArray = np.array(list(NonRedundantSet))
print(NonRedundantArray)

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
print('lengths of PairwiseSimilaritiesDict ',  len(PairwiseSimilaritiesDict.keys()))
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
SimilarityDataFrame = pd.DataFrame(data=ScoreMatrix,    # values
            index=NonRedundantSet,
            columns=NonRedundantSet)
#print(ScoreMatrix)
#print(SimilarityDataFrame)
print(datetime.datetime.now())
print('ScoreMatrix created')

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
# print(dendrogram(SimilarityMatrixArray))
# print(to_tree(SimilarityMatrixArray))
FclusterArray = fcluster(SimilarityMatrixArray, 0.2, criterion='distance', depth=2, R=None, monocrit=None)
#FclusterArray = fclusterdata(ScoreMatrix, 0.2, criterion = 'distance')
#print(FclusterArray)


# dendro = dendrogram(
#     SimilarityMatrixArray,
#     labels = NonRedundantArray,
#     #truncate_mode='lastp',  # show only the last p merged clusters
#     #p=30,  # show only the last p merged clusters
#     leaf_rotation=90.,
#     leaf_font_size=12.,
#     show_contracted=True,  # to get a distribution impression in truncated branches
# )

# plt.title('Hierarchical Clustering Dendrogram')
# plt.xlabel('sample index or (cluster size)')
# plt.ylabel('distance')
# dendrogram(
#     SimilarityMatrixArray,
#     labels = NonRedundantArray,
#     #truncate_mode='lastp',  # show only the last p merged clusters
#     #p=30,  # show only the last p merged clusters
#     leaf_rotation=90.,
#     leaf_font_size=12.,
#     show_contracted=True,  # to get a distribution impression in truncated branches
# )
# plt.show()
#
# plot_tree(dendro)

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