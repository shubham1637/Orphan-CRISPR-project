import re
from collections import Counter

GeneralTableD1FileName = "/panfs/pan1/orphancrispr/GeneralIslandsClustersTable_d1.csv"
GeneralTableD1txtFileName = "/panfs/pan1/orphancrispr/GeneralIslandsClustersTable_d1.txt"
GeneralTableD05FileName = "/panfs/pan1/orphancrispr/GeneralIslandsClustersTable_d0.5.csv"
GeneralTableD05txtFileName = "/panfs/pan1/orphancrispr/GeneralIslandsClustersTable_d0.5.txt"
CountedGeneFreqsWithinCluster05FileName = "/panfs/pan1/orphancrispr/CountedGeneFreq_withinCluster_d0.5.txt"
CountedGeneFreqsWithinCluster1FileName = "/panfs/pan1/orphancrispr/CountedGeneFreq_withinCluster_d1.txt"
SortedCountedGeneFreqsWithinCluster05FileName = "/panfs/pan1/orphancrispr/CountedGeneFreq_sorted_withinCluster_d0.5.txt"
SortedCountedGeneFreqsWithinCluster1FileName = "/panfs/pan1/orphancrispr/CountedGeneFreq_sorted_withinCluster_d1.txt"
SortedCountedGeneFreqsWithinCluster1WithSizesFileName  = "/panfs/pan1/orphancrispr/CountedGeneFreq_sorted_withinCluster_d1_.txt"
SortedCountedGeneFreqsWithinCluster05WithSizesFileName  = "/panfs/pan1/orphancrispr/CountedGeneFreq_sorted_withinCluster_d0.5_.txt"

ClustersWithGenes = {}
ClustersWithSizes = {}
count = 0
for line in open(GeneralTableD05FileName, "r"):
    count += 1
    if count == 1:
        continue
    ClusterID = line[:-1].split(',')[2].replace('"', '')
    ClustersWithSizes[ClusterID] = [line[:-1].split(',')[10], line[:-1].split(',')[12]]
    if ClusterID not in ClustersWithGenes:
        ClustersWithGenes[ClusterID] = line[:-1].split(',')[9].replace('"','').split(';') # string with gene identificators
    else:
        ClustersWithGenes[ClusterID] += line[:-1].split(',')[9].replace('"','').split(';')

print(ClustersWithSizes)

DictWithFreq = {}
for key in ClustersWithGenes:
    DictWithFreq[int(key.replace('"', ''))] = Counter(ClustersWithGenes[key])

with open(CountedGeneFreqsWithinCluster05FileName, "w") as File:
    for i in range(1,10000):
        if i in DictWithFreq:

            File.write(str(i) + '\t')
            for gene in DictWithFreq[i]:
                File.write(gene + ':' + str(DictWithFreq[i][gene]) + ' ')
            File.write('\n')


with open(SortedCountedGeneFreqsWithinCluster1FileName, "w") as File:
    for line in open(CountedGeneFreqsWithinCluster05FileName, "r"):
        GenesWithFreq = {}
        File.write(str(line[:-1].split('\t')[0]) + '\t')
        for gene in line[:-1].split('\t')[1].split(' '):
            if gene != '':
                GenesWithFreq[int(gene.split(':')[1])] = gene.split(':')[0]
        for key in sorted(GenesWithFreq.keys(),reverse=True):
            File.write(GenesWithFreq[key] + ':' + str(key) + ' ')
        File.write('\n')


with open(SortedCountedGeneFreqsWithinCluster05WithSizesFileName, "w") as File:
    for line in open(SortedCountedGeneFreqsWithinCluster05FileName, "r"):
        File.write(str(line[:-1].split('\t')[0]) + '\t' + str(ClustersWithSizes[str(line[:-1].split('\t')[0])][0]).replace('"', '') + '\t'
                   + str(ClustersWithSizes[str(line[:-1].split('\t')[0])][1]).replace('"', '')  + '\t' + line[:-1].split('\t')[1] + '\n')