import datetime
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, average, fcluster, fclusterdata, to_tree

NonRedundantSequencesFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Nonredundant_repeats.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName80 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_0.8.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName75 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_0.75.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName70 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_0.7.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName65 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_0.65.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName60 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_0.6.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName55 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_0.55.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName50 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_0.5.txt"



NonRedundantSet = set()
for line in open(NonRedundantSequencesFileName, "r"):
    NonRedundantSet.add(line[:-1])

print(datetime.datetime.now(), ' Open the DistanceMatrix from a file')
DistanceMatrixUpperDiag = np.loadtxt('/panfs/pan1/orphancrispr/ClusteringRepeats/0703/DistanceMatrixUpperDiag.txt',delimiter = ',')
#
# print(datetime.datetime.now(), ' Start of hierarchical clustering with distance = 0.2')
# SimilarityMatrixArray = linkage(DistanceMatrixUpperDiag, method = 'single')
# FclusterArray = fcluster(SimilarityMatrixArray, 0.2, criterion='distance', depth=2, R=None, monocrit=None)
# print(datetime.datetime.now(), ' End of hierarchical clustering with distance = 0.2')
#
# with open(NonRedundantRepeatsIDwithClusterLabelsFileName80, "w") as File:
#     i=0
#     for name in NonRedundantSet:
#        File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
#        i += 1

# print(datetime.datetime.now(), ' Start of hierarchical clustering with distance = 0.25')
# SimilarityMatrixArray = linkage(DistanceMatrixUpperDiag, method = 'single')
# FclusterArray = fcluster(SimilarityMatrixArray, 0.25, criterion='distance', depth=2, R=None, monocrit=None)
# print(datetime.datetime.now(), ' End of hierarchical clustering with distance = 0.25')
#
# with open(NonRedundantRepeatsIDwithClusterLabelsFileName75, "w") as File:
#     i=0
#     for name in NonRedundantSet:
#        File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
#        i += 1

# print(datetime.datetime.now(), ' Start of hierarchical clustering with distance = 0.3')
# SimilarityMatrixArray = linkage(DistanceMatrixUpperDiag, method = 'single')
# FclusterArray = fcluster(SimilarityMatrixArray, 0.3, criterion='distance', depth=2, R=None, monocrit=None)
# print(datetime.datetime.now(), ' End of hierarchical clustering with distance = 0.3')
#
# with open(NonRedundantRepeatsIDwithClusterLabelsFileName70, "w") as File:
#     i=0
#     for name in NonRedundantSet:
#        File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
#        i += 1

# print(datetime.datetime.now(), ' Start of hierarchical clustering with distance = 0.35')
# SimilarityMatrixArray = linkage(DistanceMatrixUpperDiag, method = 'single')
# FclusterArray = fcluster(SimilarityMatrixArray, 0.35, criterion='distance', depth=2, R=None, monocrit=None)
# print(datetime.datetime.now(), ' End of hierarchical clustering with distance = 0.35')
#
# with open(NonRedundantRepeatsIDwithClusterLabelsFileName65, "w") as File:
#     i=0
#     for name in NonRedundantSet:
#        File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
#        i += 1

# print(datetime.datetime.now(), ' Start of hierarchical clustering with distance = 0.4')
# SimilarityMatrixArray = linkage(DistanceMatrixUpperDiag, method = 'single')
# FclusterArray = fcluster(SimilarityMatrixArray, 0.4, criterion='distance', depth=2, R=None, monocrit=None)
# print(datetime.datetime.now(), ' End of hierarchical clustering with distance = 0.4')
#
# with open(NonRedundantRepeatsIDwithClusterLabelsFileName60, "w") as File:
#     i=0
#     for name in NonRedundantSet:
#        File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
#        i += 1

# print(datetime.datetime.now(), ' Start of hierarchical clustering with distance = 0.45')
# SimilarityMatrixArray = linkage(DistanceMatrixUpperDiag, method = 'single')
# FclusterArray = fcluster(SimilarityMatrixArray, 0.45, criterion='distance', depth=2, R=None, monocrit=None)
# print(datetime.datetime.now(), ' End of hierarchical clustering with distance = 0.45')
#
# with open(NonRedundantRepeatsIDwithClusterLabelsFileName55, "w") as File:
#     i=0
#     for name in NonRedundantSet:
#        File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
#        i += 1

print(datetime.datetime.now(), ' Start of hierarchical clustering with distance = 0.5')
SimilarityMatrixArray = linkage(DistanceMatrixUpperDiag, method = 'single')
FclusterArray = fcluster(SimilarityMatrixArray, 0.5, criterion='distance', depth=2, R=None, monocrit=None)
print(datetime.datetime.now(), ' End of hierarchical clustering with distance = 0.5')

with open(NonRedundantRepeatsIDwithClusterLabelsFileName50, "w") as File:
    i=0
    for name in NonRedundantSet:
       File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
       i += 1