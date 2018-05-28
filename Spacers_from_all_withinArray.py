import re
import itertools
import io
import subprocess

SpacersBaseFileName =  "/panfs/pan1.be-md.ncbi.nlm.nih.gov/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/Spacers.fna"
OrphanArraysFileName = "/home/utkinai2/Project1/orphans_with_ORF_full-part.txt"
Orphan99ArraysFileName = "/home/utkinai2/Project1/99newOrphans.txt"
IdentifiedArraysFileName = "/home/utkinai2/Project1/identified_with_ORF_full-part.txt"
SpacersFastafromFullOrphansFileName = "/home/utkinai2/Project1/SpacersFasta_fromFullOrphans.txt"
SpacersFastafromPartOrphansFileName = "/home/utkinai2/Project1/SpacersFasta_fromPartOrphans.txt"
SpacersFastafromIdentifiedFileName = "/home/utkinai2/Project1/SpacersFasta_fromIdentified_2.txt"
SpacersFrom99OrphansFileName = "/home/utkinai2/Project1/99newOrphans_spacers.txt"
BlastSpacersFrom99OrphansFileName = "/home/utkinai2/Project1/99newOrphans_spacers.hits"

BlastSpacersFileName = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/orphancrispr/BlastResult/FullOrphan_Spacers.hits"

ResultingTableFileName = "/home/utkinai2/Project1/Spacers_afterblast_from_99newOrphans.txt"
#ClustersWithMajorTypeFileName = "/home/utkinai2/Project1/SpacersFasta_fromClusters/Clusters_with_0.85-1proportion.txt"
#SpacersFastaFileName = "/home/utkinai2/Project1/SpacersFasta_fromClusters/tmp.txt"

# creating dictionary with keys = clusters, values =  and string containing ID of bacteria, start and end coordinates of CRISPR-array
# in a format corresponding to the Spacers.fna file (ID_start_end_..), so that all spacers for certain array could be easily extracted
Bacteria= {'Orphans': [], 'Identified': []}
for line in open(Orphan99ArraysFileName, "r"):
    LineValues = line[:-1].split(' ')
    Bacteria['Orphans'] += [[str(LineValues[2].replace('"', '') + '_' + LineValues[3] + '_' + LineValues[4]), LineValues[5]]]
# for line in open(IdentifiedArraysFileName, "r"):
#     LineValues = line[:-1].split(' ')
#     Bacteria['Identified'] += [[str(LineValues[0].replace('"', '') + '_' + LineValues[1] + '_' + LineValues[2]), LineValues[3]]]


# this part is hidden because running takes much time and necessary file was already created
# BUT SHOULD BE RAN ONE TIME
# with open(SpacersFastafromFullOrphansFileName, "w") as SpacersFastafromFullOrphansFile:
#     with open(SpacersFastafromPartOrphansFileName, "w") as SpacersFastafromPartOrphansFile:
#with open(SpacersFastafromIdentifiedFileName, "w") as SpacersFastafromIdentifiedFile:
#             for orphan in Bacteria['Orphans']:
#                 with open(SpacersBaseFileName, "r") as f:
#                     for line1, line2 in itertools.zip_longest(*[f] * 2):
#                         Match = re.findall(str(orphan[0]), line1)
#                         if Match:
#                             print(Match)
#                             if orphan[1] == "Full":
#                                 SpacersFastafromFullOrphansFile.write(line1 + line2)
#                             if orphan[1] == "Partial":
#                                 SpacersFastafromPartOrphansFile.write(line1 + line2)
#     for identified in Bacteria['Identified']:
#         with open(SpacersBaseFileName, "r") as f:
#             for line1, line2 in itertools.zip_longest(*[f] * 2):
#                 Match = re.findall(str(identified[0]), line1)
#                 if Match:
#                     print(Match)
#                     SpacersFastafromIdentifiedFile.write(line1 + line2)




#then run blastn on all three files
# subprocess.call("blastn  -query /home/utkinai2/Project1/SpacersFasta_fromClusters/Spacers_Cluster_" + str(number) + ".txt -subject /home/utkinai2/Project1/SpacersFasta_fromClusters/Spacers_Cluster_" +
#                     str(number) + ".txt" + " -word_size 8 -dust no -task blastn -outfmt '7 qseqid sseqid slen sstart send evalue qseq sseq qstart qend bitscore score staxids sscinames stitle'"
#                     " > BlastSpacers_Cluster_" + str(number) + ".hits",
#                     shell=True)

# analysis of hits, search of identical spacers within one array
IDsSpacerNumber = {}
for line in open(SpacersFrom99OrphansFileName, "r"):
    if line[0] == ">":
        LineValues = line[:-1].split('_')
        IDsSpacerNumber[str(LineValues[0].strip('>') + '_' + str(LineValues[1])+ '_' + str(LineValues[2]))] = [LineValues[3]]


Array = []
for line in open(BlastSpacersFrom99OrphansFileName, "r"):
    if line[0] == '#':
        if line[2:7] == "Query":
            LineValues = line[:-1].split(' ')
            LineVal = LineValues[2].split('_')
            Array = [str(LineVal[0] + '_' + LineVal[1] + '_' + LineVal[2]), LineVal[3]]
        else:
            continue
    else:
        if Array != [] and re.match(Array[0], line[:-1].split('\t')[1]):
            if line[:-1].split('\t')[1].split('_')[3] != Array[1]:
                IDsSpacerNumber[Array[0]] += [Array[1]]

with open(ResultingTableFileName, "w") as ResultingTableFile:
    for key in IDsSpacerNumber.keys():
        IDsSpacerNumber[key] = [int(IDsSpacerNumber[key][0]), len(list(set(IDsSpacerNumber[key][1:])))]
        print( IDsSpacerNumber[key])
        if IDsSpacerNumber[key][1] == 0 or IDsSpacerNumber[key][1] == IDsSpacerNumber[key][0]:
            ResultingTableFile.write(key.split('_')[0] + '\t' + key.split('_')[1] + '\t' + key.split('_')[2] + '\t' +
                                 str(IDsSpacerNumber[key][0]) + '\t' + str(IDsSpacerNumber[key][1]) + '\n')
        else:
            ResultingTableFile.write(key.split('_')[0] + '\t' + key.split('_')[1] + '\t' + key.split('_')[2] + '\t' +
                                     str(IDsSpacerNumber[key][0]) + '\t' + str(int(IDsSpacerNumber[key][1] + 1)) + '\n')   #because 1 matching spacer = 2 spacers in array are equal.
            #that is why +1 to number of identical spacers, so it make sense in terms of whole array
#print(IDsSpacerNumber)
print(len(IDsSpacerNumber.keys()))