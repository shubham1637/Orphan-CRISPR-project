import numpy as np
import pandas as pd

CRISPRsFileName = "/home/utkinai2/Project1/lefties.csv"
FileWithIDandContigs = "/home/utkinai2/Project1/all1603.pp.txt"
# CRISPRsFileName = "C:/Users/utkinai2/Desktop/Ira/identified.csv"
# FileWithIDandContigs = "C:/Users/utkinai2/Desktop/Ira/all1603.pp.txt"
LeftCRISPRsWithinContigFileName = "/home/utkinai2/Project1/leftovers_filter.txt"
count = 0
Offset = 2000

ContigIDSizes = {}
for line in open(FileWithIDandContigs, "r"):
    LineValues = line[:-1].split("\t")
    if line[0] == "#":
        continue
    ContigIDSizes[LineValues[1]] = int(LineValues[4])


with open(LeftCRISPRsWithinContigFileName, "w") as LeftCRISPRsWithinContigFile:
    for line in open(CRISPRsFileName, "r"):
        count += 1
        if count < 2:
            continue  # header is skipped now
        LineValues = line[:-1].split(",")
        if LineValues[2].strip("\"") in ContigIDSizes:
            if (int(LineValues[3]) < Offset) or (int(LineValues[4]) > ContigIDSizes[LineValues[2].strip("\"")] - Offset):
                continue
            LeftCRISPRsWithinContigFile.write(line)




#     count1 = 0
#     LineValues = line[:-1].split(" ")
#     for row in open(CRISPRsFileName, "r"):
#         count1 +=1
#         if count1 <2:
#             continue
#         RowValues = row[:-1].split(",")
#         if RowValues[2] == LineValues[1]:

# DataContigs = pd.read_csv(FileWithIDandContigs, sep='\t', low_memory= False, names = ["V1", "V2", "V3", "V4", "V5", "V6", "V7"])
# DataContigs = DataContigs.iloc[1:, :]
# # CotigEnd = DataContigs.iloc[:, "V5"]
# CRISPRsFile = pd.read_csv(CRISPRsFileName, sep = ',',low_memory= False)# names = ["V1", "V2", "V3", "V4", "V5", "V6", "V7","V8", "V9", "V10", "V4", "V5", "V6", "V7"] )
# #print(CRISPRsFile.tail())
# count = 0
# # print(CRISPRsFile.loc[104, 1], ' ', int(DataContigs.loc[int(CRISPRsFile[104, 1]), 4]))
# for i in range(1, CRISPRsFile.shape[0]):
#     if (CRISPRsFile.loc[i, "V2"] in DataContigs.loc[2:, "V2"]):
#     #     # print(CRISPRsFile.loc[i, "V2"], ' ', int(DataContigs.loc[CRISPRsFile[i, "V2"], "V5"]))
#     #     if (int(CRISPRsFile.loc[i, "V3"]) < 2000) | (int(CRISPRsFile.loc[i, "V4"]) >  int(DataContigs.loc[CRISPRsFile[i, "V2"], "V5"]) - 2000):
#         count += 1
#     print(count)

# print(CRISPRsFile.loc[1, "V2"])