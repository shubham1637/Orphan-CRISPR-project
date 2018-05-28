import re

SourceFileName = "/home/utkinai2/Project1/merged_all.csv"
ClusteredFileName = "/home/utkinai2/Project1/ClusteringRepeats/BlastClust/all_blastclust_0.96_out.txt"
MappedWithClusterFileName = "/home/utkinai2/Project1/ClusteringRepeats/BlastClust/all_blastclust_clusters_96%.txt"

CRISPRarrayInfoDict = {}
cluster = 0
count = 0
for Line in open(SourceFileName, "r"):
    count += 1
    if count < 2:
        continue
    LineVal = Line[:-1].split(",")
    CRISPRarrayInfoDict[int(LineVal[0].strip('"'))] = LineVal[2] + "\t" + LineVal[3] + "\t" + LineVal[4] + "\t" +\
                                                 LineVal[6] + "\t" + LineVal[10] + "\t" + LineVal[11] + "\t" + LineVal[15] # LineVal[15] for orphan + identified only


with open(MappedWithClusterFileName, "w") as File:
    for line in open(ClusteredFileName, "r"):
        cluster += 1
        LineValues = line[:-2].split(' ')
        for SeqNumber in LineValues:
            if int(SeqNumber) in CRISPRarrayInfoDict:
                print(cluster, ' yeah!')
                File.write(str(cluster) + '\t' + CRISPRarrayInfoDict[int(SeqNumber)] + '\n')

