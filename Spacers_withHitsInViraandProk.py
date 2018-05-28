import itertools

SpacerSourceFileName = '/home/utkinai2/Project1/SpacersFasta_from_Orphans.txt'
FilteredOrphansFileName = "/home/utkinai2/Project1/Orphans_finalfiltered.txt"
SpacersFilteredSourceFileName = '/panfs/pan1/orphancrispr/SpacersFromArraysCluster/Spacers_from_Orphans.txt'
ViralHitsFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/Orphans_vs_viruses.hits"
ViralHitsFilteredFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/FinalOrphans_vs_viruses_filtered95-95.hits"
ArraysWithViralHitsFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/ArraysWithViralHits_95-95.txt"
OrphansWithProkHitsFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/Orphans_filtered_SpacerHits.hits"
FinalOrphansWithProkHitsFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/FinalOrphans_filtered_SpacerHits.hits"

OrphanID = set()
for line in open(FilteredOrphansFileName, "r"):
    OrphanID.add(str( line[:-1].split(' ')[2].replace('"','') + '_' +line[:-1].split(' ')[3] + '_' + line[:-1].split(' ')[4] ))
# Set = set()
# with open(SpacersFilteredSourceFileName, "w") as File:
#     with open(SpacerSourceFileName, "r") as f:
#         for line1, line2 in itertools.zip_longest(*[f] * 2):
#             ID = line1[1:-1].split('_')[0] + '_' + line1[1:-1].split('_')[1] + '_' + line1[1:-1].split('_')[2]
#             if ID in OrphanID:
#                 Set.add(ID)
#                 File.write(line1 + line2)
# print(len(Set))




# SpacerIDandLength = {}
# with open(SpacersFilteredSourceFileName, "r") as f:
#     for line1, line2 in itertools.zip_longest(*[f] * 2):
#         SpacerIDandLength[line1[1:-1]] = len(line2[:-1])

# ArrayIDWithViralHit = {}
# with open(ViralHitsFilteredFileName, "w") as File:
#     for line in open(ViralHitsFileName, "r"):
#         if line[0] == "#":
#             continue
#
#         LineValues = line[:-1].split('\t')
#         ID = LineValues[0].split('_')[0] + '_' + LineValues[0].split('_')[1] + '_' + LineValues[0].split('_')[2] + '_'\
#              + LineValues[0].split('_')[3] + '_' + LineValues[0].split('_')[4]
#         if ID not in SpacerIDandLength:
#             continue
#
#         seq1, seq2 = LineValues[6], LineValues[7]
#         Identity = 0
#         for i in range(len(min(seq1, seq2))):
#             if seq1[i] == seq2[i]:
#                 Identity += 1/len(min(seq1, seq2))
#         Coverage =  len(seq1) / SpacerIDandLength[ID]
#         if Coverage >= 0.95 and Identity >= 0.95:
#             print(Coverage, ' ', Identity, ' ', ID)
#             ArrayIDWithViralHit[LineValues[0]] = [LineValues[14]]
#             File.write(line)
#
# with open(ArraysWithViralHitsFileName, "w") as file:
#     for key in ArrayIDWithViralHit:
#         file.write(key + '\t' + ' '.join(ArrayIDWithViralHit[key]) + '\n')


with open(FinalOrphansWithProkHitsFileName, "w") as File:
    for line1 in open(OrphansWithProkHitsFileName, "r"):
        ID = line1[:-1].split('_')[0] + '_' + line1[:-1].split('_')[1] + '_' + line1[:-1].split('_')[2]
        if ID in OrphanID:
            File.write(line1)