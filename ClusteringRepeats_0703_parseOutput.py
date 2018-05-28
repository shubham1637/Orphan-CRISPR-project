import io
import re

AllArraysInfoFileName = "/home/utkinai2/Project1/merged_all.csv"

CRISPRinfoByID = {}
for line in open(AllArraysInfoFileName,"r"):
    LineValues = line[:-1].split(",")
    CRISPRinfoByID[str(LineValues[2] + '_' + LineValues[3]+ '_' + LineValues[4])] = str(LineValues[5] + '\t' + LineValues[6] + '\t' +
                                                                                        LineValues[10] + '\t' + LineValues[11] + '\t' +
                                                                                        LineValues[15])

ClustersDict = {}
Similarity = [0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]
for threshold in Similarity:
    with io.open('/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Clusters_ofNonRedundantRepeats_' + str(threshold) + '.txt', "w") as File:
        print(threshold)
        for line in io.open("/panfs/pan1/orphancrispr/ClusteringRepeats/0703/NonredundantRepeats_with_ClusterLabels_" + str(threshold) + ".txt" , "r"):
            LineValues  = line[:-1].split('\t')
            if LineValues[0] in CRISPRinfoByID:
                File.write(str(LineValues[1]) + '\t' + LineValues[0] + '\t' + CRISPRinfoByID[LineValues[0]] + '\n')

        #     if LineValues[1] in ClustersDict:
        #         ClustersDict[LineValues[1]] += [LineValues[0]]
        #     else:
        #         ClustersDict[LineValues[1]] = [LineValues[0]]
        #
        # for cluster in sorted(ClustersDict):
        #     File.write('Cluster ' + str(cluster) + '\n' + str(ClustersDict[cluster]).replace('"', '') + '\n')