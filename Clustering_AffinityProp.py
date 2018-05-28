from sklearn.cluster import AffinityPropagation
import numpy as np
import pandas as pd
from sklearn import metrics
# from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, average, fcluster, fclusterdata
import numpy as np
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist

RepeatsSimilaritiesFileName = "/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/AllRepeatsVsAllRepeats_fix.hits"
IDwithLabelsFileName = "/panfs/pan1/orphancrispr/IDwithClusterLabels.txt"
NameListFile = "/panfs/pan1/orphancrispr/NameList2.txt"
AllFileName = "/home/utkinai2/Project1/merged_all_with_clusters.txt"
PairwiseSimilaritiesDictFileName = "/home/utkinai2/Project1/Pairwise_similaritiesDict_without_2.txt"
AssignedDictFileName = "/home/utkinai2/Project1/AssignedIDs_with_Score2.txt"
# NameList0 = []
# for line in open(AllFileName, "r"):
#     LineValues = line[:-1].split('\t')
#     NameList0.append(str(LineValues[2] +'_' + LineValues[3] +'_' + LineValues[4]))
# print(len(NameList0))

# SimilarityDataFrame = pd.DataFrame(index=NameList0, columns=NameList0)
# SimilarityDataFrame = SimilarityDataFrame.fillna(0)

NameList = []
AssignedDict = {}
count = 0
PairwiseSimilaritiesDict = {}
# NonRedundantSet = set()
for line in open(RepeatsSimilaritiesFileName, "r"):
    count += 1
    LineValues = line[:-1].split('\t')
    Name1 = LineValues[0]
    Name2 = LineValues[1]
    Score = int(LineValues[3])/(2*int(LineValues[2]))
    # print(Score)
    if Score == 1:
        # print(count, ' YES')
        if Name1 not in AssignedDict.keys():
            AssignedDict[Name1] = Name2
        else:
            AssignedDict[Name1] += Name2
    if Name1 not in AssignedDict.values():
        NameList.append(Name1)
    # SimilarityDataFrame.loc[Name1, Name2] = Score
    # if count%500 == 0:
    #     print(count)
    else:
        PairwiseSimilaritiesDict[Name1, Name2] = Score

with open(PairwiseSimilaritiesDictFileName,"w" ) as f:
    for key in PairwiseSimilaritiesDict:
       f.write(key + '\t' + PairwiseSimilaritiesDict[key] + '\n')

with open(AssignedDictFileName, "w") as AssignedFile:
    for key in AssignedDict:
        AssignedFile.write(key + '\t' + AssignedDict[key] + '\n')

#
print(len(NameList))
print(len(PairwiseSimilaritiesDict.keys()))

with open(NameListFile,"w") as NameFile:
    for name in NameList:
        NameFile.write(name + '\t')
#

# SimilarityDataFrame = pd.DataFrame(index=NonRedundantSet, columns=NonRedundantSet)
# SimilarityDataFrame = SimilarityDataFrame.fillna(0)
#
# # for line in open(RepeatsSimilaritiesFileName, "r"):
# #     count += 1
# #     LineValues = line[:-1].split('\t')
# #     Name1 = LineValues[0]
# #     Name2 = LineValues[1]
# #     Score = int(LineValues[3])/(2*int(LineValues[2]))
# #     SimilarityDataFrame.loc[Name1, Name2] = Score
#
#
#
#
# for name1 in NonRedundantSet:
#     for name2 in NonRedundantSet:
#         if (name1, name2) in PairwiseSimilaritiesDict:
#             SimilarityDataFrame.loc[name1,name2] = PairwiseSimilaritiesDict[name1,name2]
#
# # 1 print("SimilarityDataFrame completed")
# # print(SimilarityDataFrame.iloc[1:10,1:10])
#
# # compute the similarity matrix
# # af = AffinityPropagation( affinity='precomputed', verbose=True)
# # aff = AffinityPropagation( affinity='precomputed', verbose=True).fit(SimilarityDataFrame)
# # af.fit_predict(SimilarityDataFrame)
# # print (af.labels_, af.cluster_centers_indices_)
# #
# # cluster_centers_indices = af.cluster_centers_indices_
# # labels = af.labels_
# # print(len(labels))
# #
# # n_clusters_ = len(cluster_centers_indices)
# #
# # print('Estimated number of clusters: %d' % n_clusters_)
# # NonRedundantSet = list(NonRedundantSet)
# # labels = list(labels)
# # with open(RepresentativeIDwithLabelsFileName, "w") as File:
# #         File.write(NonRedundantSet + '\t' + labels  + '\n')
#
#
# Z = linkage(SimilarityDataFrame, method = 'average')
# FclusterArray = fcluster(Z, 0.95, criterion='inconsistent', depth=2, R=None, monocrit=None)
# print(FclusterArray)
# c, coph_dists = cophenet(Z, pdist(SimilarityDataFrame))
# print(c)
# with open(RepresentativeIDwithLabelsFileName, "w") as File:
#     i=0
#     for name in NonRedundantSet:
#        File.write(name + '\t' + FclusterArray[i] + '\n')
#        i += 1

# print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
# print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
# print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
# print("Adjusted Rand Index: %0.3f"
#       % metrics.adjusted_rand_score(labels_true, labels))
# print("Adjusted Mutual Information: %0.3f"
#       % metrics.adjusted_mutual_info_score(labels_true, labels))
# print("Silhouette Coefficient: %0.3f"
#       % metrics.silhouette_score(X, labels, metric='sqeuclidean'))
# # Plot result
# import matplotlib.pyplot as plt
# from itertools import cycle
#
# plt.close('all')
# plt.figure(1)
# plt.clf()
#
# colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
# for k, col in zip(range(n_clusters_), colors):
#     class_members = af.labels_ == k
#     cluster_center = SimilarityDataFrame[af.cluster_centers_indices_[k]]
#     plt.plot(SimilarityDataFrame[class_members, 0], SimilarityDataFrame[class_members, 1], col + '.')
#     plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
#              markeredgecolor='k', markersize=14)
#     for x in SimilarityDataFrame[class_members]:
#         plt.plot([cluster_center[0], x[0]], [cluster_center[1], x[1]], col)
#
# plt.title('Estimated number of clusters: %d' % n_clusters_)
# plt.show()