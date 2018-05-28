import datetime
import numpy as np
import itertools
from scipy.cluster.hierarchy import dendrogram, linkage, average, fcluster, fclusterdata, to_tree

NonRedundantSequencesFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Nonredundant_repeats.txt"
RepeatsAssignedFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/RepeatsIdentical.txt"
RepeatsSimilaritiesFileName = "/panfs/pan1/orphancrispr/AllRepeatsVsAllRepeats_fix.hits"
AllArraysInfoFileName ="/home/utkinai2/Project1/merged_all.csv"
# IDwithLabelsFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/All_SS_withClusterLabels_forRepresentatives_0.9.txt"
# AllIDwithLabelsFileName = "/home/utkinai2/Project1/ClusteringRepeats/All_SS_withClusterLabels_all_0.9.txt"

RepeatLengthsByID = {}
for line in open(AllArraysInfoFileName,"r"):
    LineValues = line[:-1].split(",")
    RepeatLengthsByID[str(LineValues[2] + '_' + LineValues[3] + '_' + LineValues[4])] = len(LineValues[6])

def ReverseComplement(seq):
   complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
   reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
   return reverse_complement

print(datetime.datetime.now())

SequenceToClustDict = dict()
ID = ""
for Line in open("/panfs/pan1/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/AllRepeats.fna"):
   if Line[0] == ">":
       ID = Line[1:-1]
   else:
       Repeat = Line[:-1]
       ReverseComplementRepeat = ReverseComplement(Repeat)

       if Repeat in SequenceToClustDict:
           SequenceToClustDict[Repeat].append(ID)
       elif ReverseComplementRepeat in SequenceToClustDict:
           SequenceToClustDict[ReverseComplementRepeat].append(ID)
       else:
           SequenceToClustDict[Repeat] = [ID]

#creating
with open(RepeatsAssignedFileName, "w") as NonRedundantFile:
    for Repeat in SequenceToClustDict:
        NonRedundantFile.write(Repeat + '\t' + ' '.join(SequenceToClustDict[Repeat]) + '\n')

NonRedundantSet = set()
for line in open(RepeatsAssignedFileName,"r"):
    RepresentativeID = line[:-1].split('\t')[1].split(' ')[0]
    NonRedundantSet.add(RepresentativeID)

with open(NonRedundantSequencesFileName ,"w") as f:
    for case in NonRedundantSet:
        f.write(case + '\n')


