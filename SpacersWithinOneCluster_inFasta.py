import re
import itertools
import io
import subprocess

SpacersBaseFileName =  "/panfs/pan1.be-md.ncbi.nlm.nih.gov/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/Spacers.fna"
ClustersFileName = "/home/utkinai2/Project1/all_with_0.92similar_clusters.txt"
SpacersFastafromClustersFileName = "/home/utkinai2/Project1/SpacersFasta_fromClusters/SpacersFasta_fromClusters.txt"
ClustersWithMajorTypeFileName = "/home/utkinai2/Project1/SpacersFasta_fromClusters/Clusters_with_0.85-1proportion.txt"
#SpacersFastaFileName = "/home/utkinai2/Project1/SpacersFasta_fromClusters/tmp.txt"

# creating dictionary with keys = clusters, values =  and string containing ID of bacteria, start and end coordinates of CRISPR-array
# in a format corresponding to the Spacers.fna file (ID_start_end_..), so that all spacers for certain array could be easily extracted
BacteriaDict1 = {}
for line in open(ClustersFileName, "r"):
    LineValues = line[:-1].split('\t')
    if LineValues[0] not in BacteriaDict1.keys():
        Bacteria = []  # array for saving unique bacteria within this cluster
        BacteriaDict1[LineValues[0]] = [str(LineValues[1] + '_' + LineValues[2] + '_' + LineValues[3])]
        Bacteria.append(LineValues[6])
    else:
        if LineValues[6] not in Bacteria:
            Bacteria.append(LineValues[6])
            BacteriaDict1[LineValues[0]] += [str(LineValues[1] + '_' + LineValues[2] + '_' + LineValues[3])]
        else:
            continue

# creating a file containing fasta of spacers of all arrays within each cluster, separated by line "Cluster <number_of_cluster>"
# part below is hidden because running this part takes a lot of time and necessary file is already created

# with open(SpacersFastafromClustersFileName, "w") as SpacersFastafromClustersFile:
#     for cluster in BacteriaDict1.keys():
#         SpacersFastafromClustersFile.write('Cluster' + '\t' + str(cluster) + '\n')
#         for array in BacteriaDict1[cluster]:
#             # SpacersFastafromClustersFile.write('test')
#             with open(SpacersBaseFileName, "r") as f:
#                 for line1, line2 in itertools.zip_longest(*[f] * 2):
#                     Match = re.findall(str(array.strip('[]')), line1)
#                     if Match:
#                         print(Match)
#                         SpacersFastafromClustersFile.write(line1 + line2)
ClustersList = []
# creating a list with clusters, where the proportion of any major CEISPR type to all types present within cluster is more than 0.85
for line in open(ClustersWithMajorTypeFileName, "r"):
    LineValues = line[:-1].split(' ')
    ClustersList.append(int(LineValues[1]))
print(ClustersList)

# creating files containing fasta of all spacers within one cluster
# THIS IS WORKING PART
# hidden because running takes a lot of time; has been run once, all necessary files created
flag = 0
count = 0
# for number in ClustersList:
#     with io.open("Spacers_Cluster_" + str(number) + ".txt" , "w") as TempFile:
#         for line in open(SpacersFastafromClustersFileName, "r"):
#             LineValues = line[:-1].split('\t')
#             if flag == 1:
#                 TempFile.write(line + '\n')
#             if LineValues[0] == "Cluster" and int(LineValues[1]) == number:
#                 flag = 1
#                 #print(flag)
#             if LineValues[0] == "Cluster" and int(LineValues[1]) != number:
#                 flag = 0
                #print(flag)


for number in ClustersList:
    print(number)
    subprocess.call("blastn  -query /home/utkinai2/Project1/SpacersFasta_fromClusters/Spacers_Cluster_" + str(number) + ".txt -subject /home/utkinai2/Project1/SpacersFasta_fromClusters/Spacers_Cluster_" +
                    str(number) + ".txt" + " -word_size 8 -dust no -task blastn -outfmt '7 qseqid sseqid slen sstart send evalue qseq sseq qstart qend bitscore score staxids sscinames stitle'"
                    " > BlastSpacers_Cluster_" + str(number) + ".hits",
                    shell=True)

