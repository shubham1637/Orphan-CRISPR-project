
NonRedundantRepeatsWithNumbersFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Nonredundant_repeats_with_numbers.txt"
AllArraysInfoFileName ="/home/utkinai2/Project1/merged_all.csv"
HclustLabelsFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/Hclust/ClusterLabels_1.txt"
NonRedundantSequencesFileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/0703/Nonredundant_repeats_with_numbers.txt"
# ClusteredThreshold100FileName = "/panfs/pan1/orphancrispr/ClusteringRepeats/Hclust/Threshold100.txt"


CRISPRinfoByID = {}
for line in open(AllArraysInfoFileName,"r"):
    LineValues = line[:-1].split(",")
    CRISPRinfoByID[str(LineValues[2] + '_' + LineValues[3]+ '_' + LineValues[4])] = str(LineValues[5] + '\t' + LineValues[6] + '\t' +
                                                                                        LineValues[10] + '\t' + LineValues[11] + '\t' +
                                                                                        LineValues[15])
NonRedundantSet = []
for line in open(NonRedundantSequencesFileName, "r"):
    NonRedundantSet.append(line[:-1].split('\t')[0])

count = -2
Threshold = [100, 95, 90, 85,84,83,82,81, 80, 75, 70, 65, 60, 55, 50,45,40,35,30,25,20,15,10,5,0]
for line in open(HclustLabelsFileName, "r"):
    count += 1
    if count == -1:
        continue
    with open("/panfs/pan1/orphancrispr/ClusteringRepeats/Hclust/Threshold_" + str(Threshold[count]) + "_1.txt","w") as File:
        LineValues = line[:-1].split(' ')
        i = 0
        print( str(NonRedundantSet[i]).split('_')[0] + '\t' + str(NonRedundantSet[i]).split('_')[1] + '\t' +
                       str(NonRedundantSet[i]).split('_')[2] + '\t' + CRISPRinfoByID[str(NonRedundantSet[i])])
        for cluster in LineValues[1:]:
            File.write(str(cluster) + '\t' + str(NonRedundantSet[i]).split('_')[0] + '\t' + str(NonRedundantSet[i]).split('_')[1] + '\t' +
                       str(NonRedundantSet[i]).split('_')[2] + '\t' + CRISPRinfoByID[str(NonRedundantSet[i])] + '\n')
            i += 1

