import itertools
import re

#
# Labels = ["merge", "labels", "height"]
#
# for uno in Labels:
#     with open("//panfs/pan1/orphancrispr/IslandsCluster/wg_output_" + str(uno) + ".txt","w") as File:
#         count = 0
#         print(uno)
#         for line in open("//panfs/pan1/orphancrispr/IslandsCluster/wg_" + str(uno) + ".txt","r"):
#             count +=1
#             if count == 1:
#                 continue
#             LineValues = line[:-1].split(' ')
#             File.write('\t'.join(LineValues[1:]).replace('"', '') + '\n')

# sum = 0
# count = 0
# wg = []
# for line in open('//panfs/pan1/orphancrispr/IslandsCluster/weights2003.output', "r"):
#     count +=1
#     sum += float(line[:-1].split('\t')[1])
#     wg.append(line[:-1].split('\t')[0])
# print(sum, ' ', count)



# with open('/panfs/pan1/orphancrispr/Spacers_from_finalOrphans.txt', "w") as File:
#     for line in open('/panfs/pan1/orphancrispr/OrphanArrays_with_Spacers.txt', "r"):
#         count = 0
#         ID = line[:-1].split('\t')[0]
#         Spacers = line[:-1].split('\t')[1].split(',')
#         for spacer in Spacers:
#             count+=1
#             File.write('>' + ID + '_' + str(count) + '_' + "spacer" + '\n' +
#                        spacer + '\n')

#
# with open("/panfs/pan1/orphancrispr/SpacersFromArraysCluster/IdentifiedSpacerHits_filtered90-90.hits", "w") as FilteredSpacerHitsFile:
#     for line in open("/panfs/pan1/orphancrispr/SpacersFromArraysCluster/Spacers_from_Identified_AllagainstAll.hits" ,"r"):
#         if line[0] == "#":
#             continue
#         LineValues = line[:-1].split('\t')
#         if float(LineValues[6]) >= 90.0 and float(LineValues[7]) >= 90.0:
#             #print('YES' , ID[0], ' ', identity, '% ', coverage, '%')
#             FilteredSpacerHitsFile.write(line)
# # Set1 = set()
# Set2 = set()
# List = []
# for line in open("//panfs/pan1/orphancrispr/IslandsCluster/IslandsID_spacers_withNumber.txt", "r"):
#     Set1.add(line[:-1].split('\t')[0] )
# for line in open("//panfs/pan1/orphancrispr/IslandsCluster/IslandsID_pfam_withNumber.txt", "r"):
#     List.append(line[:-1].split('\t')[0])
#     Set2.add(line[:-1].split('\t')[0])
# List = list(Set2 - Set1)
# print(List)
# #
# SpacersBaseFileName =  "/panfs/pan1.be-md.ncbi.nlm.nih.gov/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/Spacers.fna"
#
# with open('//panfs/pan1/orphancrispr/IslandsCluster/Missing31_arraySpacers.txt',"w") as File:
#     with open(SpacersBaseFileName, "r") as f:
#         for line1, line2 in itertools.zip_longest(*[f] * 2):
#             for ID in List:
#                 Match = re.findall(ID, line1)
#                 if Match:
#                     print(Match)
#                     File.write(line1 + line2)