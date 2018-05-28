import re

SourceFileName = "/home/utkinai2/Project1/merged_all.csv"
ClusteredFileName = "/home/utkinai2/Project1/ClusteringRepeats/CDHIT_0103/All_repeats_90%.sorted"
MappedWithClusterFileName = "/home/utkinai2/Project1/ClusteringRepeats/CDHIT_0103/All_with_0.9similar_clusters.txt"
# NewWithClusterFileName = "/home/utkinai2/Project1/All_with_1.0similar_clusters_.txt"
#ReverseRepeatsFileName = "/home/utkinai2/Project1/reverse_repeats.txt"

SeqNumber = 0
Percentage = 0
Cluster = 0
Number = 0
with open(MappedWithClusterFileName, "w") as MappedWithClusterFile:
#with open(ReverseRepeatsFileName, "w") as ReverseFile:
    for Line in open(ClusteredFileName, "r"):
        if Line[0] == ">":
            if re.search(r'[0-9]{1,}',Line[:-1].split(" ")[1]):
                Cluster = re.search(r'[0-9]{1,}',Line[:-1].split(" ")[1]).group(0)
                continue
        LineValues = Line[:-1].split(" ")
        if re.search(r'[0-9]{1,}', LineValues[1]):
            SeqNumber = re.search(r'[0-9]{1,}', LineValues[1]).group(0)
        for line in open(SourceFileName,"r"):
            LineValues2 = line[:-1].split(",")
            if re.search(r'[0-9]{1,}', LineValues2[0]):
                Number = re.search(r'[0-9]{1,}', LineValues2[0]).group(0)
            if Number == SeqNumber:
                print(Cluster)
              #  ReverseFile.write(">" + LineValues2[2].replace('"', ''))
                MappedWithClusterFile.write(str(Cluster) + '\t' + LineValues2[2].replace('"', '') + '\t' +
                                            str(LineValues2[3]) + '\t' + str(LineValues2[4]) + '\t' +
                                            LineValues2[6].replace('"', '') + '\t' + LineValues2[10].replace('"', '') + '\t' +
                                            LineValues2[11].replace('"', '') + '\t' + LineValues2[15].replace('"', '') + '\n')                                     # + LineValues2[15].replace('"', '')  - for merged_id_orphans file
#                break
# with open(NewWithClusterFileName, "w") as NewFile:
#     for line in open(MappedWithClusterFileName, "r"):
#         line = line.replace('Leftovers', '\tOrphans')
#         line = line.replace('Identified', '\tIdentified')
#         NewFile.write(line)