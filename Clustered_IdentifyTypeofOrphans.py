import re
from operator import itemgetter
import collections
from collections import Counter
ClustersFilteredFileName = "/home/utkinai2/Project1/identified_with_1.0similar_clusters.txt"
TypesInClustersFileName = "/home/utkinai2/Project1/identified_types_in_1.0similar_clusters.txt"

Clusters = {}
n = 0
Bacteria = str()
PreviousCluster = -1
#with open(TypesInClustersFileName, "w") as TypesInClustersFile:
for line in open(ClustersFilteredFileName, "r"):
    LineValues = line[:-1].split("\t")
    Bacteria = LineValues[6].replace('"', '').split("_")[0]
    #Class = re.search(r'[a-zA-Z]{6,10}', LineValues[9]).group(0)
    # if Class == "Identified" and LineValues[7].replace('"','') == "Unidentified":
    #     continue
    if LineValues[5] == "Unidentified":
        continue
    if int(LineValues[0]) != PreviousCluster:
        Clusters[int(LineValues[0])] = [LineValues[5].replace('"', ''), Bacteria]
                      # Class]]
    if int(LineValues[0]) == PreviousCluster:
        Clusters[int(LineValues[0])] += [LineValues[5].replace('"', ''), Bacteria]
                                       # Class]]
    PreviousCluster = int(LineValues[0])

#print(Clusters)
#print(Clusters[21])
counts = []

# for cluster in Clusters.keys():
#     counter = collections.Counter(Clusters[cluster])
#     print(counter)

print(collections.Counter(Clusters[2]))

#print(Clusters[21])
#
ClusterComponentsCounted = []
ClusterComponents = {}

# for i in Clusters.keys():
#     Clusters[i] = sorted(Clusters[i], key=itemgetter(0))
#     Lastcase = Clusters[i][0]
#     ClusterComponentsCounted = []
#     Bact = []
#     count = 0
#     for case in Clusters[i]:
#         if case[0] == Lastcase[0]:
#             count += 1
#             # if case[1] not in Bact:
#             #     Bact.append(case[1])
#         else:
#             ClusterComponentsCounted.append([Lastcase, count])
#             count = 1
#             Lastcase = case
#     ClusterComponentsCounted.append([Lastcase, count])
#     ClusterComponents[i] = ClusterComponentsCounted
#
# with open(TypesInClustersFileName, "w") as TypesInClustersFile:
#     for key in ClusterComponents.keys():
#         TypesInClustersFile.write('Cluster' + '\t' + str(key) +'\n')
#         for value in ClusterComponents[key]:
#             for elem in value:
#                 TypesInClustersFile.write(str(elem) + '\n')



