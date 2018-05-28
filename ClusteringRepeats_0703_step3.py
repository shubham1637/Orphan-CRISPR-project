import datetime
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, average, fcluster, fclusterdata, to_tree
import pandas as pd

FilteredRepeatsSimilaritiesFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/AllRepeatsVsAllRepeats_fix_filtered_Dist.hits"
NonRedundantSequencesFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Nonredundant_repeats.txt"
NonRedundantRepeatsWithNumbersFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Nonredundant_repeats_with_numbers.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName100 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_100.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName95 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_95.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName90 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_90.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName85 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_85.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName80 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_80.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName75 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_75.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName70 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_70.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName65 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_65.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName60 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_60.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName55 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_55.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName50 = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_50.txt"


TestFilteredRepeatsSimilaritiesFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Test_filtered.hits "
# tri = np.array([[0,0.1,0.3,0.05, 1.0],[1.0,0,1.0,0.1, 1.0],[1.0,0.16,0,0.5, 0.43], [0.3,1.0,1.0,0, 0.17], [0.3,0.4,1.0,0, 1.0]])
# tri_upper_diag = np.triu(tri, k=0)
# print(tri_upper_diag)
# np.savetxt('//panfs/pan1/orphancrispr/ClusteringRepeats/0703/test_nparray.txt', tri_upper_diag, delimiter = ',')
# a = np.loadtxt('//panfs/pan1/orphancrispr/ClusteringRepeats/0703/test_nparray.txt',delimiter = ',')
# print(a)
# Z = linkage(a, method="single")
# FclusterArray = fcluster(Z, 0.2, criterion='distance', depth=2, R=None, monocrit=None)
# print(FclusterArray)

NonRedundantSet = set()
for line in open(NonRedundantSequencesFileName, "r"):
    NonRedundantSet.add(line[:-1])

print(datetime.datetime.now())
count = 0
NonRedundantSetWithNumbers = {}
with open(NonRedundantRepeatsWithNumbersFileName, "w") as NonRedundantRepeatsWithNumbersFile:
    for Name in NonRedundantSet:
        NonRedundantSetWithNumbers[Name] = count
        NonRedundantRepeatsWithNumbersFile.write(Name + '\t' + str(count) +'\n')
        count += 1
print(datetime.datetime.now())

count1=0
DistanceMatrix = np.ones((len(NonRedundantSet), len(NonRedundantSet)))
#np.full((2, 2), 10)

for line in open(FilteredRepeatsSimilaritiesFileName, "r"):
    count1 += 1
    LineValues = line[:-1].split('\t')
    Name1 = LineValues[0]
    Name2 = LineValues[1]
    DistanceMatrix[NonRedundantSetWithNumbers[Name1],NonRedundantSetWithNumbers[Name2]] = LineValues[2]

print(datetime.datetime.now(), ' Distance Matrix is created')

#



df = pd.DataFrame(DistanceMatrix)
df.to_csv("/panfs/pan1/orphancrispr/ClusteringRepeats/DistanceMatrix_1303.csv")
#DistanceMatrix = DistanceMatrix.flatten()
# print(DistanceMatrix)
#DistanceMatrixUpperDiag = np.triu(DistanceMatrix, k=0)
#DistanceMatrix = np.reshape(DistanceMatrix, (330512400,1))
#print(DistanceMatrix[30:50,])

# DistanceMatrixUpperDiag = np.triu(DistanceMatrix, k=0)
# print( DistanceMatrixUpperDiag[100:130, 18090:18120])

# RUN FOLLOWING ROW ONLY ONCE CAUSE IT TAKES TIME!
# np.savetxt('/panfs/pan1/orphancrispr/ClusteringRepeats/0703/DistanceMatrix_1.txt', DistanceMatrix, delimiter = ',')

# print(datetime.datetime.now(), ' Start of hierarchical clustering distance = 0')
# SimilarityMatrixArray = linkage(DistanceMatrix, method = 'single')
# FclusterArray = fcluster(SimilarityMatrixArray, 0.0, criterion='distance', depth=2, R=None, monocrit=None)
# print(datetime.datetime.now(), ' End of hierarchical clustering')
#
# with open(NonRedundantRepeatsIDwithClusterLabelsFileName100, "w") as File:
#     i=0
#     for name in NonRedundantSet:
#        File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
#        i += 1

# print(datetime.datetime.now(), ' Start of hierarchical clustering with distance = 0.05')
# SimilarityMatrixArray = linkage(DistanceMatrix, method = 'single')
# FclusterArray = fcluster(SimilarityMatrixArray, 0.05, criterion='distance', depth=2, R=None, monocrit=None)
# print(datetime.datetime.now(), ' End of hierarchical clustering with distance = 0.05')
#
# with open(NonRedundantRepeatsIDwithClusterLabelsFileName95, "w") as File:
#     i=0
#     for name in NonRedundantSet:
#        File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
#        i += 1

# print(datetime.datetime.now(), ' Start of hierarchical clustering with distance = 0.1')
# SimilarityMatrixArray = linkage(DistanceMatrix, method = 'single')
# FclusterArray = fcluster(SimilarityMatrixArray, 0.1, criterion='distance', depth=2, R=None, monocrit=None)
# print(datetime.datetime.now(), ' End of hierarchical clustering with distance = 0.1')
#
# with open(NonRedundantRepeatsIDwithClusterLabelsFileName90, "w") as File:
#     i=0
#     for name in NonRedundantSet:
#        File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
#        i += 1
#
