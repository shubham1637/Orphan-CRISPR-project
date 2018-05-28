import re
import itertools

# SpacerIDWithNumberFileName = "/panfs/pan1/orphancrispr/SpacersFromArraysCluster/SpacerID_withNumber_2003.txt"
SpacerIDWithNumberFileName = "/panfs/pan1/orphancrispr/SpacersFromArraysCluster/IdentifiedSpacerID_withNumber.txt"
SpacerSourceFileName = '/panfs/pan1/orphancrispr/SpacersFromArraysCluster/Spacers_from_Identified.txt'
AssignedToNonRedundantFileName = "/panfs/pan1/orphancrispr/SpacersFromArraysCluster/NonRedundantSpacers_withAssigned_Identified_wg.txt"
SpacersWithClusterAssignedFileName = "/panfs/pan1/orphancrispr/SpacersFromArraysCluster/Spacer_with_ClusterNumber_h0.5_Identified.txt"
# ArrayWithSpacersFileName =  "/panfs/pan1/orphancrispr/SpacersFromArraysCluster/OrphanArrays_with_Spacers.txt"
ArrayWithClusteredSpacersFileName =  "/panfs/pan1/orphancrispr/SpacersFromArraysCluster/IdentifiedArrays_with_ClusteredSpacers.txt"


SpacerInfoByID = {}
with open(SpacerSourceFileName,"r") as f:
    for line1, line2 in itertools.zip_longest(* [f] * 2):
        ID = line1[1:-1]
        SpacerInfoByID[ID] = line2[:-1]

AssignedSpacerIDs = {}
for line in open(AssignedToNonRedundantFileName, "r"):
    NonRedundantSpacer = line[:-1].split('\t')[0]
    for assigned in line[:-1].split('\t')[1].replace('[', '').replace(']','').split(','):
        if NonRedundantSpacer in AssignedSpacerIDs:
            AssignedSpacerIDs[NonRedundantSpacer] += [assigned.strip('"').strip("'").strip(" '")]
        else:
            AssignedSpacerIDs[NonRedundantSpacer] = [assigned.strip('"').strip("'").strip(" '")]


Threshold = [0, 0.5, 0.1] #, 0.2, 0.3, 0.4 ,0.5, 0.6, 0.7, 0.8, 0.9, 1]
for threshold in Threshold:
    with open("/panfs/pan1/orphancrispr/SpacersFromArraysCluster/Output_HclustSpacers_Identified_h" + str(threshold) + ".txt","w") as File:
        count = 0
        for line in open("/panfs/pan1/orphancrispr/SpacersFromArraysCluster/HclustSpacers_Identified_h" + str(threshold) + ".txt" ,"r"):
            #print(threshold)
            count +=1
            if count == 1:
                continue

            LineValues = line[:-1].split(' ')
            ID = LineValues[0].strip('"')
            cluster = LineValues[1]
            File.write(str(cluster) + '\t' + ID + '\t' + SpacerInfoByID[ID] + '\n')

NonRedundantSpacersWithCluster = {}

for line in open("/panfs/pan1/orphancrispr/SpacersFromArraysCluster/Output_HclustSpacers_Identified_h0.5.txt" , "r"):
    NonRedundantSpacersWithCluster[line[:-1].split('\t')[1]] = line[:-1].split('\t')[0]

with open(SpacersWithClusterAssignedFileName, "w") as File:
    for spacer in AssignedSpacerIDs:
        Cluster = 'Cluster_' + str(NonRedundantSpacersWithCluster[spacer])
        File.write(spacer + '\t' + Cluster + '\n')
        for assignedspacer in AssignedSpacerIDs[spacer]:
            File.write(assignedspacer + '\t' + Cluster + '\n')

ArrayIDwithClusteredSpacers = {}
for line in open(SpacersWithClusterAssignedFileName, "r"):
    ArrayID = line[:-1].split('\t')[0].split('_')[0] + '_' + line[:-1].split('\t')[0].split('_')[1] + '_' + line[:-1].split('\t')[0].split('_')[2]
    ClusterOfSpacer = line[:-1].split('\t')[1]
    if ArrayID in ArrayIDwithClusteredSpacers:
        if ClusterOfSpacer not in ArrayIDwithClusteredSpacers[ArrayID]:
            ArrayIDwithClusteredSpacers[ArrayID] += [ClusterOfSpacer]
    else:
        ArrayIDwithClusteredSpacers[ArrayID] = [ClusterOfSpacer]
print(len(ArrayIDwithClusteredSpacers))


with open(ArrayWithClusteredSpacersFileName, "w") as File:
    for array in ArrayIDwithClusteredSpacers:
        File.write(array + '\t' + ','.join(ArrayIDwithClusteredSpacers[array]) + '\n' )
