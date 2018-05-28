

IslandsAnnotationFileName = "/panfs/pan1/orphancrispr/IslandsCluster/Islands_pfam.ann_clust.xls"
GeneralTableD1FileName = "/panfs/pan1/orphancrispr/GeneralIslandsClustersTable_d1.txt"
# GeneralTableD1ProteinsFileName = "/panfs/pan1/orphancrispr/GeneralIslandsClustersTable_d1_.txt"
LabeledIslandsAnnotationFileName = "/panfs/pan1/orphancrispr/IslandsCluster/Islands_pfam_labeled.ann_clust"
NewLabeledIslandsAnnotationFileName = "/panfs/pan1/orphancrispr/IslandsCluster/Islands_pfam_clusterlabeled.ann_clust"
CountedGeneFreqsWithinCluster1FileName = "/panfs/pan1/orphancrispr/CountedGeneFreq_sorted_withinCluster_d1.txt"
CountedGeneFreqsWithinCluster1ProtFileName = "/panfs/pan1/orphancrispr/CountedGeneFreq_sorted_withinCluster_d1_prot.txt"

PfamProtein = {}
IDandCluster = {}

for line in open(GeneralTableD1FileName, "r"):
    IDandCluster[line[:-1].split('\t')[2][:-2] + '_' + line[:-1].split('\t')[3] + '_' + line[:-1].split('\t')[4]] = str(line[:-1].split('\t')[1])
print(len(IDandCluster))


HeaderAndCluster = {}
with open(LabeledIslandsAnnotationFileName, "w") as File:
    for line in open(IslandsAnnotationFileName, "r"):
        LineValues = line[:-1].split('\t')

        if LineValues[0] == "===":
            header = line
            File.write(line)
        # if LineValues[0] == "===":
        #     File.write(LineValues[0] + '\t' + "Neighborhood_" + str(IDandCluster[LineValues[4]]) + '\t' +
        #                LineValues[2] + '\t' + LineValues[3] + '\t' + LineValues[4] + '\t' + LineValues[5] +
        #                '\t' + LineValues[6] + '\t' + LineValues[7] + '\n')
        else:
            if LineValues[6] != "Unknown" and LineValues[6] != "CRISPR":
                count = 0
                for identifier in LineValues[6].split(','):
                    if len(LineValues[8].split(',')) <= len(LineValues[6].split(',')):
                        PfamProtein[identifier] = LineValues[8].split(',')[count]
                    count += 1

            if LineValues[6] != "CRISPR":
                File.write(line)
            if LineValues[6] == "CRISPR":
                if LineValues[4] + '_' + str(LineValues[1].split('..')[0]) + '_' + str(LineValues[1].split('..')[1]) in IDandCluster:
                    File.write(str(LineValues[0]) + '\t' + LineValues[1] + '\t' + LineValues[2] + '\t' + LineValues[3] + '\t' + LineValues[4] + '\t' + str(LineValues[5]) +
                           '\t' + LineValues[6] + '\t' + LineValues[7] + '\t' + LineValues[8] + '\t' + str(LineValues[9]) + '\t' + LineValues[10] + '\t' +
                               "Neighborhood_" + IDandCluster[LineValues[4] + '_' + str(LineValues[1].split('..')[0]) + '_' + str(LineValues[1].split('..')[1])]
                               + '\t' + LineValues[12] + '\t' + LineValues[13] + '\t' + LineValues[14] + '\n')
                    HeaderAndCluster[header] = "Neighborhood_" + IDandCluster[LineValues[4] + '_' + str(LineValues[1].split('..')[0]) + '_' + str(LineValues[1].split('..')[1])]
                else:
                    File.write(line)
lll = 0
with open(NewLabeledIslandsAnnotationFileName, "w") as File:
    for line in open(LabeledIslandsAnnotationFileName, "r"):
        if line in HeaderAndCluster:
            LineValues = line[:-1].split('\t')
            File.write(LineValues[0] + '\t' + HeaderAndCluster[line] + '\t' + LineValues[2] + '\t' + LineValues[3] + '\t' +
                       LineValues[4] + '\t' + LineValues[5] + '\t' + LineValues[6] + '\t' + LineValues[7] + '\n')
        else:
            File.write(line)

        if line[0] == '=' and line not in HeaderAndCluster:
            lll += 1
print(lll)

with open(CountedGeneFreqsWithinCluster1ProtFileName, "w") as File:
    for line in open(CountedGeneFreqsWithinCluster1FileName, "r"):
        LineValues = line[:-1].split('\t')
        File.write(LineValues[0] + '\t' + LineValues[1] + '\t' + LineValues[2] + '\t')
        for prot in LineValues[3].split(' '):
            if prot.split(':')[0] in PfamProtein:
                File.write(prot.split(':')[0] + '(' + PfamProtein[prot.split(':')[0]] + '):' + prot.split(':')[1] + ' ')
            else:
                File.write(prot + ' ')
        File.write('\n')
