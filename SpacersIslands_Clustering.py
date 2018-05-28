import re
import itertools

SpacersBaseFileName =  "/home/utkinai2/Project1/Spacers_from_Orphans.txt"
IslandIDwithCRISPRNeighborsFixedFileName = "/panfs/pan1/orphancrispr/IslandsCluster/Islands_pfam_3NeighborsSet_fixedID.txt"
#IslandIDwithCRISPRNeighborsFixedFileName = "/panfs/pan1/orphancrispr/IslandsCluster/h0.2_2cluster_pfam.txt"
ArrayWithSpacersFileName =  "/panfs/pan1/orphancrispr/OrphanArrays_with_Spacers.txt"
#ArrayWithClusteredSpacersFileName =  "/panfs/pan1/orphancrispr/SpacersFromArraysCluster/OrphanArrays_with_ClusteredSpacers.txt"
#ArrayWithSpacersFileName =  "/panfs/pan1/orphancrispr/h0.2_2cluster_spacers.txt"


# PART FOR CLUSTERING BY SPACERS (SEQUENCES)

def ReverseComplement(seq):
   complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
   reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
   return reverse_complement

SpacersFromArray = {}

#run this once
for line in open(IslandIDwithCRISPRNeighborsFixedFileName , "r"):
    ID = line[:-1].split(' ')[2].replace('"','') + '_' + line[:-1].split(' ')[3] + '_' + line[:-1].split(' ')[4]
    print(ID)
    with open(SpacersBaseFileName, "r") as f:
        for line1, line2 in itertools.zip_longest(*[f] * 2):
            Match = re.findall(ID, line1)
            if Match:
                # print(Match)
                if ID in SpacersFromArray:
                    SpacersFromArray[ID] += [line2[:-1], ReverseComplement(line2[:-1])]
                else:
                    SpacersFromArray[ID] = [line2[:-1], ReverseComplement(line2[:-1])]

print('Length of SpacersFromaArray dict = ', len(SpacersFromArray))
with open(ArrayWithSpacersFileName, "w") as File:
    for ID in SpacersFromArray:
        File.write(ID + '\t' + ','.join(SpacersFromArray[ID]) + '\n')


