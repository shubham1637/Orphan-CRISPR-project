import datetime
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage, average, fcluster, fclusterdata
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist

NonRedundantSequencesFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/Nonredundant_repeats.txt"
AssignedWithScore1FileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/IDsAssigned_with_Score-1.txt"
RepeatsSimilaritiesFileName = "/panfs/pan1/orphancrispr/AllRepeatsVsAllRepeats_fix.hits"
AllArraysInfoFileName ="/home/utkinai2/Project1/merged_all.csv"
IDwithLabelsFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/All_SS_withClusterLabels_forRepresentatives_0.9.txt"
AllIDwithLabelsFileName = "/home/utkinai2/Project1/ClusteringRepeats/All_SS_withClusterLabels_all_0.9.txt"

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

Counter = 0
with open(AssignedWithScore1FileName, "w") as AssignedWithScore1:
    for Repeat in SequenceToClustDict:
       Counter += 1
       AssignedWithScore1.write(str(Counter) + "\t" + Repeat + "\t" + ",".join(SequenceToClustDict[Repeat]) + '\n')

DictWithAssignedSequences = {}
NonRedundantSequenceWithID = {}
NonRedundantIDs = []
for line in open(AssignedWithScore1FileName, "r"):
    LineValues = line[:-1].split('\t')
    IDs = LineValues[-1].split(',')
    DictWithAssignedSequences[IDs[0]] = [IDs[1:]]
    NonRedundantSequenceWithID[IDs[0]] = LineValues[1]
    NonRedundantIDs.append(IDs[0])
print(len(NonRedundantSequenceWithID))

with open(NonRedundantSequencesFileName, "w") as NonRedundantSequencesFile:
    for ID in NonRedundantSequenceWithID:
        NonRedundantSequencesFile.write('>' + ID + '\n' + NonRedundantSequenceWithID[ID] + '\n')

# print(datetime.datetime.now())
# SimilarityDataFrame = pd.DataFrame(index=NonRedundantIDs, columns=NonRedundantIDs)
# SimilarityDataFrame = SimilarityDataFrame.fillna(0)
# print(datetime.datetime.now())
# print('empty DataFrame of shape ', SimilarityDataFrame.shape, ' created')


PairwiseSimilaritiesDict = {}
count = 0
for line in open(RepeatsSimilaritiesFileName, "r"):
    LineValues = line[:-1].split('\t')
    Name1 = LineValues[0]
    Name2 = LineValues[1]
    if Name1 == Name2:
        PairwiseSimilaritiesDict[Name1, Name2] = 1
        continue
    Score = int(LineValues[3]) * min(RepeatLengthsByID[Name1], RepeatLengthsByID[Name2]) / (
    max(RepeatLengthsByID[Name2], RepeatLengthsByID[Name1]) * 2 * int(LineValues[2]))
    if Name1 in NonRedundantSequenceWithID and Name2 in NonRedundantSequenceWithID:
        count+=1
        print(count)


# Assigning scores from PairwiseSimilaritiesDict to the corresponding indices to the zero-matrix
ScoreMatrix = np.zeros((len(NonRedundantSet), len(NonRedundantSet)))
for (Name1, Name2) in PairwiseSimilaritiesDict:
    ScoreMatrix[NonRedundantSetInDigits[Name1], NonRedundantSetInDigits[Name2]] = PairwiseSimilaritiesDict[Name1, Name2]
    ScoreMatrix[NonRedundantSetInDigits[Name2], NonRedundantSetInDigits[Name2]] = PairwiseSimilaritiesDict[Name1, Name2]
#print(ScoreMatrix)
#print(SimilarityDataFrame)
print(datetime.datetime.now())

with open(PairwiseSimilaritiesDictFileName,"w" ) as f:
    for key in PairwiseSimilaritiesDict:
       f.write(str(key) + '\t' + str(PairwiseSimilaritiesDict[key]) + '\n')

SimilarityMatrixArray = linkage(ScoreMatrix, method = 'average')
FclusterArray = fcluster(SimilarityMatrixArray, 0.1, criterion='distance', depth=2, R=None, monocrit=None)
#print(FclusterArray)

#print(c)

with open(IDwithLabelsFileName, "w") as File:
    i=0
    for name in list(NonRedundantIDs):
       File.write(str(name) + '\t' + str(FclusterArray[i]) + '\n')
       i += 1

print(datetime.datetime.now())
print('FclusterArray created')

ClustersDict = {}
with open(AllIDwithLabelsFileName, "w") as File:
    for line in open(IDwithLabelsFileName, "r"):
        TempList = set()
        LineValues = line[:-1].split('\t')
        if LineValues[0] in NonRedundantSequenceWithID:
            for i in range(len(NonRedundantSequenceWithID[LineValues[0]])):
                TempList.add(str(NonRedundantSequenceWithID[LineValues[0]][i]))
            print(TempList)
            TempList = list(TempList)
        if len(TempList) != 0:
            if LineValues[1] in ClustersDict:
                ClustersDict[LineValues[1]] += [LineValues[0],str([str(TempList[i]) for i in range(len(TempList))])[2:-2]]
            else:
                ClustersDict[LineValues[1]] = [LineValues[0], str([str(TempList[i]) for i in range(len(TempList))])[2:-2]]
        else:
            if LineValues[1] in ClustersDict:
                ClustersDict[LineValues[1]] += [LineValues[0]]
            else:
                ClustersDict[LineValues[1]] = [LineValues[0]]

    for cluster in sorted(ClustersDict):
        File.write('Cluster ' + str(cluster) + '\n' + str(ClustersDict[cluster]).replace('"','') + '\n')

print('DONE.')
print(datetime.datetime.now())