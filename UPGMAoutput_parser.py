import re
import itertools

SourceFileName = "/home/utkinai2/Project1/merged_all.csv"
ClusteredFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/WithClusterLabels_all_0.9_dist_0203.txt"
MappedWithClusterFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/TableWithClusters_0.9similar_UPGMA.txt"

CRISPRarrayInfoDict = {}
count=0
for Line in open(SourceFileName, "r"):
    count += 1
    if count < 2:
        continue
    LineVal = Line[:-1].split(",")
    ID = str(LineVal[2] + '_' + LineVal[3] + '_' + LineVal[4]).replace('"','')
    print(ID)
    CRISPRarrayInfoDict[ID] = LineVal[6] + "\t" + LineVal[10] + "\t" + LineVal[11] + "\t" + LineVal[15] #+ '\t' + LineVal[16]# LineVal[15] for orphan + identified only and LineVal[16] for reverserepeats

with open(MappedWithClusterFileName, "w") as File:
   with open(ClusteredFileName, "r") as f:
       for line1, line2 in itertools.zip_longest(*[f] * 2):
           Cluster = line1[8:-1]
           LineValues = str(line2)[1:-2].split(',')
           for ID in LineValues:
                    ID = ID.replace(" ", '').replace("'", '')
                    print(ID)
                    if ID in CRISPRarrayInfoDict:
                        print('!')
                        File.write(str(Cluster) + '\t' + str(ID.replace('_', '\t')) + '\t' + str(CRISPRarrayInfoDict[ID]) + '\n')

