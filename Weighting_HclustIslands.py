

IDandWeightFileName = "/panfs/pan1/orphancrispr/IslandsCluster/weights2003.output"

IDandWeight = {}
for line in open(IDandWeightFileName, "r"):
    IDandWeight[line[:-1].split('\t')[0]] = line[:-1].split('\t')[1]


Heights = [0, 0.2, 0.5,0.75, 1,1.25, 1.5, 1.6, 1.65, 1.7, 1.75 , 1.8, 1.9, 2, 2.25, 2.5, 3, 3.5, 4, 4.5, 5]
for uno in Heights:
    with open("/panfs/pan1/orphancrispr/IslandsCluster/HclustIslands" + str(uno) + "_weighted_2003.txt", "w") as File:
        count = 0
        print(uno)
        for line in open("/panfs/pan1/orphancrispr/IslandsCluster/HclustIslands" + str(uno) , "r"):
            count +=1
            if count == 1:
                continue

            LineValues = line[:-1].split(' ')
            ID = LineValues[2].replace('"', '') + '_' + LineValues[3] + '_' + LineValues[4]
            if ID in IDandWeight:
                 File.write(line[:-1].replace('"','') + ' ' + str(IDandWeight[ID]) + '\n')

