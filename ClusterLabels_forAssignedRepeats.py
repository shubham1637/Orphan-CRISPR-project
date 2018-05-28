

RepeatsAssignedFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/RepeatsIdentical.txt"
NonRedundantRepeatsIDwithClusterLabelsFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/Hclust/Threshold_85_1.txt"
AllRepeatsWithClusterLabelsFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/AllRepeats_with_ClusterLabels_0.85.txt"
ArraysWithViralHitsFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/ArraysWithViralHits_95-95.txt"
GeneralTableFileName = "/panfs/pan1/orphancrispr/GeneralIslandsClustersTable_d1.txt"
TableFileName = "/panfs/pan1/orphancrispr/IslandsCluster/IslandsClusters_weightedBySpacers_dist1.txt"

NonredRepeatsWithlabels = {}
for line in open(NonRedundantRepeatsIDwithClusterLabelsFileName, "r"):
    NonredRepeatsWithlabels[line[:-1].split('\t')[1] + '_' + line[:-1].split('\t')[2] + '_' + line[:-1].split('\t')[3]]  = line[:-1].split('\t')[0]

with open(AllRepeatsWithClusterLabelsFileName, "w") as file:
    for line in open(RepeatsAssignedFileName, "r"):
        for ID in line[:-1].split('\t')[1].split(' '):
            file.write(ID + '\t' + str(NonredRepeatsWithlabels[line[:-1].split('\t')[1].split(' ')[0]]) + '\n')

AllRepeatsWithLabels = {}
for line in open(AllRepeatsWithClusterLabelsFileName, "r"):
    AllRepeatsWithLabels[line[:-1].split('\t')[0]] = line[:-1].split('\t')[1]

IDandViralHit = {}
for line in open(ArraysWithViralHitsFileName, "r"):
    ID = line[:-1].split('\t')[0].split('_')[0] + '_' + line[:-1].split('\t')[0].split('_')[1] + '_' + line[:-1].split('\t')[0].split('_')[2]
    IDandViralHit[ID] = line[:-1].split('\t')[1]
print(IDandViralHit)

count =0
with open(GeneralTableFileName, "w") as Table:
    for line in open(TableFileName, "r"):
        LineValues = line[:-1].replace('"','').split(' ')
        count +=1
        if count == 1:
            Table.write("RepeatCluster" + '\t' + LineValues[1] + '\t' + LineValues[2] +'\t' + LineValues[3] +'\t' + LineValues[4] +
                        '\t' + LineValues[5] + '\t' + LineValues[6] + '\t' + LineValues[7] + '\t' + LineValues[8] + '\t' + LineValues[9] + '\t' + LineValues[10] +
                        '\t' + LineValues[11] +  '\t' + "ViralMatch" + '\n')
        if count > 1:
            ID = LineValues[3] + '_' + LineValues[4] + '_' + LineValues[5]
            if ID in IDandViralHit:
                Table.write( AllRepeatsWithLabels[ID] +'\t' + LineValues[2] + '\t' + LineValues[3] +'\t' + LineValues[4] +
                        '\t' + LineValues[5] + '\t' + LineValues[6] + '\t' + LineValues[7] + '\t' + LineValues[8] + '\t' + LineValues[9] + '\t' + LineValues[10] +
                        '\t' + LineValues[11] + '\t' + LineValues[12] + '\t' + IDandViralHit[ID] + '\n')
            else:
                Table.write(
                    AllRepeatsWithLabels[ID] + '\t' + LineValues[2] + '\t' + LineValues[3] + '\t' + LineValues[4] +
                    '\t' + LineValues[5] + '\t' + LineValues[6] + '\t' + LineValues[7] + '\t' + LineValues[8] + '\t' +
                    LineValues[9] + '\t' + LineValues[10] +
                    '\t' + LineValues[11] + '\t' + LineValues[12] + '\t' + ' ' + '\n')