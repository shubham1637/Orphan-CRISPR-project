import re

EcoliSpacerHitsFileName = "/home/utkinai2/Project1/EcoliRNA_first.hits"
EcoliORFsFileName = "/home/utkinai2/Project1/Ecoli_chromosome.pty"
EcoliFastaJoinedFileName = "/home/utkinai2/Project1/Ecoli_chromosome_joined.fna"
MatchesFileName = "/home/utkinai2/Project1/PAMSpacer_matches_out_1st.txt"
SpacerFilterFileName =  "/home/utkinai2/Project1/PAMSpacer_filtered.txt"
OutORFplusFileName = "/home/utkinai2/Project1/OutORFplus.txt"
OutORFminusFileName = "/home/utkinai2/Project1/OutORFminus.txt"
count=0

SpacersLocation =  {'+': [], '-': []}
# with open (SpacerFilterFileName, "w") as SpacerFilterFile:
#     for Line in open(EcoliSpacerHitsFileName, "r"):
#         if Line[0] == "#":
#             continue
#         count += 1
#         if count < 2:
#             continue  # header is skipped now
#         LineValues = Line[:-1].split('\t') # the last one is perenos stroki, thats'why -1
#         # if re.search(r'#', LineValues[0]):
#         #     continue
#         if float(LineValues[5]) < 1.03e-08 and (int(LineValues[9]) == 33) and (int(LineValues[8]) == 1) :
#             SpacerFilterFile.write(Line)



for Line in open(SpacerFilterFileName, "r"):
    LineValues = Line[:-1].split('\t')
    if int(LineValues[3]) < int(LineValues[4]):
        SpacersLocation['+'].append([int(LineValues[3]),int(LineValues[4])])
    if int(LineValues[3]) > int(LineValues[4]):
        SpacersLocation['-'].append([int(LineValues[4]),int(LineValues[3])])

print(len(SpacersLocation['-']), ' ',len(SpacersLocation['+']))

ORFs = {'+': [], '-': []}
for Line in open(EcoliORFsFileName, "r"):
    if Line[0] == "#":
        continue
    LineValues = Line[:-1].split("\t")  # the last one is perenos stroki, thats'why -1
    RangeMatch = re.search(r'([0-9]{1,}..[0-9]{1,})', LineValues[1])
    Range = RangeMatch.group(0).split("..")
    if LineValues[2] == "+":
        ORFs['+'].append([int(Range[0]), int(Range[1])])
    if LineValues[2] == "-":
        ORFs['-'].append([int(Range[0]), int(Range[1])])

print(len(ORFs['-']), ' ',len(ORFs['+']))

Matches = []
count1 = 0
count2 = 0
Counter = 0
countPlus = 0
countMinus = 0
OutORFplus = 0
OutORFminus = 0
OverlapPlus = 0
OverlapMinus = 0
for Line in open(EcoliFastaJoinedFileName,"r"):
    #print(Line[0+57: 3+57])
    for spacer in SpacersLocation['+']:
        countPlus += 1
        for i in range(len(ORFs['+']) - 1):
            if (spacer[0] > ORFs['+'][i][1]) and (spacer[1] < ORFs['+'][i+1][0]):
                OutORFplus += 1
                Matches.append([spacer[0],spacer[1], Line[spacer[0] - 3:spacer[0]], "outplus"])
                break
    for spacer in SpacersLocation['-']:
        countMinus += 1
        for i in range(len(ORFs['-']) - 1):
            if (spacer[0] > ORFs['-'][i][1]) and (spacer[1] < ORFs['-'][i+1][0]):
                OutORFminus += 1
                Matches.append([spacer[0], spacer[1], Line[spacer[1] -1 :spacer[1] -1  + 3], "outminus"])
                break
# with open (OutORFplusFileName, "w") as O:
#     O.write(OutORFplus)
# with open (OutORFminusFileName, "w") as O:
#     O.write(OutORFminus)
# print(countPlus, ' ', countMinus)
# print(count1, ' ' , count2 )
# print(OverlapPlus, ' ', OverlapMinus)
print(OutORFplus,' ', OutORFminus)
with open(MatchesFileName, "w") as MatchesFile:
     for elem in Matches:
        MatchesFile.write(str(elem[0]) + ' ' + str(elem[1]) + ' ' + str(elem[2]) + ' ' +str(elem[3]) + "\n")