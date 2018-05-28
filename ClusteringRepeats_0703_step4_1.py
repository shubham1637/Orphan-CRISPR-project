import datetime
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, average, fcluster, fclusterdata, to_tree

NonRedundantSequencesFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Nonredundant_repeats.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName90 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_0.9.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName85 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_0.85.txt"


NonRedundantSet = set()
for line in open(NonRedundantSequencesFileName, "r"):
    NonRedundantSet.add(line[:-1])

print(datetime.datetime.now(), ' Open the DistanceMatrix from a file')
DistanceMatrixUpperDiag = np.loadtxt('/panfs/pan1/orphancrispr/ClusteringRepeats/0703/DistanceMatrixUpperDiag.txt',
                                     delimiter=',')

#
# print(datetime.datetime.now(), ' Start of hierarchical clustering with distance = 0.1')
# SimilarityMatrixArray = linkage(DistanceMatrixUpperDiag, method='single')
# FclusterArray = fcluster(SimilarityMatrixArray, 0.1, criterion='distance', depth=2, R=None, monocrit=None)
# print(datetime.datetime.now(), ' End of hierarchical clustering with distance = 0.1')
#
# with open(NonRedundantRepeatsIDwithClusterLabelsFileName90, "w") as File:
#     i = 0
#     for name in NonRedundantSet:
#         File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
#         i += 1

print(datetime.datetime.now(), ' Start of hierarchical clustering with distance = 0.15')
SimilarityMatrixArray = linkage(DistanceMatrixUpperDiag, method='single')
FclusterArray = fcluster(SimilarityMatrixArray, 0.15, criterion='distance', depth=2, R=None, monocrit=None)
print(datetime.datetime.now(), ' End of hierarchical clustering with distance = 0.15')

with open(NonRedundantRepeatsIDwithClusterLabelsFileName85, "w") as File:
    i = 0
    for name in NonRedundantSet:
        File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
        i += 1
