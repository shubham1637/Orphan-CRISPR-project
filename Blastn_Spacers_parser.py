import re

EcoliSpacerHitsFileName = "/home/utkinai2/Project1/Ecoli_AS27second_naive.hits"
EcoliORFsFileName = "/home/utkinai2/Project1/Ecoli_chromosome.pty"
EcoliFastaJoinedFileName = "/home/utkinai2/Project1/Ecoli_chromosome_joined.fna"
MatchesFileName = "/home/utkinai2/Project1/PAMSpacer_matches_AS27second_naive.txt"
SpacerFilterFileName =  "/home/utkinai2/Project1/SpacerHits_AS27second_naive_filtered.txt"
OutORFplusFileName = "/home/utkinai2/Project1/OutORFplus.txt"
OutORFminusFileName = "/home/utkinai2/Project1/OutORFminus.txt"
AAGPositionsFileName = "/home/utkinai2/Project1/AAG_positions_orf_2.txt"
count=0


# This block is running long, run it once and then comment!
SpacersLocation = {'+': [], '-': []}
with open (SpacerFilterFileName, "w") as SpacerFilterFile:
    for Line in open(EcoliSpacerHitsFileName, "r"):
        if Line[0] == "#":
            continue
        #print(Line)
        count += 1
        if count < 2:
            continue  # header is skipped now
        LineValues = Line[:-1].split('\t') # the last one is perenos stroki, thats'why -1
        #print(LineValues[1], 'd', LineValues[5])
        # if re.search(r'#', LineValues[0]):
        #     continue
        if float(LineValues[5]) < 1.03e-08 and (int(LineValues[9]) == 33) and (int(LineValues[8]) == 1) :
            SpacerFilterFile.write(Line)

Spacers = []

# for Line in open(SpacerFilterFileName, "r"):
#     LineValues = Line[:-1].split('\t')
#     if int(LineValues[3]) < int(LineValues[4]):
#         SpacersLocation['+'].append([int(LineValues[3]),int(LineValues[4])])
#         Spacers.append([int(LineValues[3]),int(LineValues[4]), LineValues[7], '+'])
#     if int(LineValues[3]) > int(LineValues[4]):
#         SpacersLocation['-'].append([int(LineValues[4]),int(LineValues[3])])
#         Spacers.append([int(LineValues[4]), int(LineValues[3]),  LineValues[7], '-'])
# print(len(SpacersLocation['-']), ' ',len(SpacersLocation['+']))
# #print(Spacers)

ORFplusLength = 0
ORFminusLength = 0
ORFs = {'+': [], '-': []}
ORFs2 = []
NonCodingLengths = []
Last = 255
CodingLength = 0
NonCodingLength = 0

Last = 320000
for Line in open(EcoliORFsFileName, "r"):
    if Line[0] == "#":
        continue
    LineValues = Line[:-1].split("\t")  # the last one is perenos stroki, thats'why -1
    RangeMatch = re.search(r'([0-9]{1,}..[0-9]{1,})', LineValues[1])
    Range = RangeMatch.group(0).split("..")


    if (int(Range[0]) >= 320000) and (int(Range[1]) <= 360000):
        NonCodingLengths.append(int(Range[0]) - Last)
        NonCodingLength += int(Range[0]) - Last
        CodingLength += int(Range[1])- int(Range[0])
        ORFs2.append([int(Range[0]), int(Range[1])])
    Last = int(Range[1])
    if LineValues[2] == "+":
        ORFs['+'].append([int(Range[0]), int(Range[1])])
        ORFplusLength += int(Range[1]) - int(Range[0])
    if LineValues[2] == "-":
        ORFs['-'].append([int(Range[0]), int(Range[1])])
        ORFminusLength += int(Range[1]) - int(Range[0])
print(NonCodingLength)
print(CodingLength)
# print(ORFminusLength,' ', ORFplusLength)
# print(len(ORFs['-']), ' ',len(ORFs['+']))

AAG_Positions = []
Matches = []
countInPlus = 0
countInMinus = 0
Counter = 0
countPlus = 0
countMinus = 0
OutORFplus = 0
OutORFminus = 0
OverlapPlus = 0
OverlapMinus = 0
OutORF = 0
InORF = 0


# for Line in open(EcoliFastaJoinedFileName,"r"):
#     print('start')
#     for match in re.finditer(r'AAG', Line):
#         # LastORFEnd = 0
#         # for orf in ORFs2:
#         #     if int(match.start()) > orf[0] and int(match.start()) < orf[1]:
#         #         AAG_Positions.append([int(match.start()), "in"])
#         #         break
#         #     if (int(match.start()) > LastORFEnd) and (int(match.start()) < orf[0]):
#         #         AAG_Positions.append([int(match.start()), "out"])
#         #         break
#         LastORFend = 0
#         for orf in ORFs['+']:
#             if int(match.start()) > orf[0] and int(match.start()) < orf[1]:
#                 AAG_Positions.append([int(match.start()), "in", "+"])
#                 break
#             if (int(match.start()) > LastORFend) and (int(match.start()) < orf[0]):
#                 AAG_Positions.append([int(match.start()), "out", "+"])
#                 break
#             LastORFend = orf[1]
#     for match in re.finditer(r'CTT', Line):
#         LastORFend = 0
#         for orf in ORFs['-']:
#             if int(match.start()) > orf[0] and int(match.start()) < orf[1]:
#                 AAG_Positions.append([int(match.start()), "in", "-"])
#                 break
#             if (int(match.start()) > LastORFend) and (int(match.start()) < orf[0]):
#                 AAG_Positions.append([int(match.start()), "out", "-"])
#                 break
#             LastORFend = orf[1]

    # for spacer in Spacers:
    #     if spacer[3] == "+":
    #         ORFend = 0
    #         for orf in ORFs['+']:
    #             if (spacer[0] > orf[0]) and (spacer[1] < orf[1]):
    #                 if (Line[spacer[0] - 3:spacer[0]] == "AAG") :
    #                     InORF += 1
    #                     elem = [spacer[0],spacer[1], spacer[2],spacer[3], Line[spacer[0] - 3:spacer[0]], "in", spacer[0]-3]
    #                     if elem not in Matches:
    #                         Matches.append(elem)
    #             if (spacer[0] > ORFend) and (spacer[1] < orf[0]):
    #                 if (Line[spacer[0] - 3:spacer[0]] == "AAG"):
    #                     OutORF += 1
    #                     elem = [spacer[0], spacer[1],spacer[2],spacer[3], Line[spacer[0] - 3:spacer[0]], "out", spacer[0] - 3]
    #                     if elem not in Matches:
    #                         Matches.append(elem)
    #             ORFend = orf[1]
    #     if spacer[3] == "-":
    #         ORFend = 0
    #         for orf in ORFs['-']:
    #             if (spacer[0] > orf[0]) and (spacer[1] < orf[1]):
    #                 if Line[spacer[1]-1:spacer[1]-1 + 3] == "CTT":
    #                     InORF += 1
    #                     elem = [spacer[0],spacer[1], spacer[2],spacer[3], Line[spacer[1]-1:spacer[1]-1 + 3], "in", spacer[1]-1]
    #                     if elem not in Matches:
    #                         Matches.append(elem)
    #             if (spacer[0] > ORFend) and (spacer[1] < orf[0]):
    #                 if Line[spacer[1] - 1:spacer[1] - 1 + 3] == "CTT":
    #                     OutORF += 1
    #                     elem = [spacer[0], spacer[1],spacer[2],spacer[3], Line[spacer[1] - 1:spacer[1] -1 +3], "out", spacer[1]-1]
    #                     if elem not in Matches:
    #                         Matches.append(elem)
    #             ORFend = orf[1]
#     #print(Line[0+57: 3+57])
#     # for spacer in SpacersLocation['+']:
#     #     countPlus += 1
#     #     for orf in ORFs['+']:
#     #         if (spacer[0] > orf[0]) and (spacer[1] < orf[1]):
#     #             #if (Line[spacerstart - 3:spacerstart] == "AAG"): #or (Line[int(spacerstart[0]-3):int(spacerstart[0])] == "CTT"):
#     #                 #+57 появляется из-за того, что с 0 по 56 знаки - это заголовок, геном начинается с 57, из-за чего сдвиг
#     #                 # spacerstart.append("AAG")
#     #                 # spacerstart.append('+')
#     #             if (Line[spacer[0] - 3:spacer[0]] == "AAG") :
#     #                 countInPlus += 1
#     #                 Matches.append([spacer[0],spacer[1], Line[spacer[0] - 3:spacer[0]], "+"])
#     #                 break
#     #         if ((spacer[0] > orf[0]) and (spacer[1] > orf[1]) and (spacer[0] < orf[1])) or (
#     #                 (spacer[0] < orf[0]) and (spacer[1] < orf[1]) and (spacer[1] > orf[0])):
#     #             if (Line[spacer[0] - 3:spacer[0]] == "AAG"):
#     #                 OverlapPlus += 1
#     #                 Matches.append([spacer[0], spacer[1], Line[spacer[0] - 3:spacer[0]], "overlapplus"])
#     #                 break
#     #
#     #
#     #     for i in range(len(ORFs['+']) - 1):
#     #         if (spacer[0] > ORFs['+'][i][1]) and (spacer[1] < ORFs['+'][i + 1][0]):
#     #             print(spacer[0],' ', ORFs['+'][i][1],' ', ORFs['+'][i + 1][0])
#     #             if (Line[spacer[0] - 3:spacer[0]] == "AAG"):
#     #                 OutORFplus += 1
#     #                 Matches.append([spacer[0], spacer[1], Line[spacer[0] - 3:spacer[0]], "outplus"])
#     #                 break
#     #         if ((spacer[0] > orf[0]) and (spacer[1] > orf[1]) and (spacer[0] < orf[1])) or ((spacer[0] < orf[0]) and (spacer[1] < orf[1]) and (spacer[1] > orf[0])):
#     #             OverlapPlus += 1
#     #             Matches.append([spacer[0], spacer[1], Line[spacer[0] - 3:spacer[0]], "overlapplus"])
#     #             break
#     #         else:
#     #             OutORFplus +=1
#     #             Matches.append([spacer[0], spacer[1], Line[spacer[0] - 3:spacer[0]], "outplus"])
#     # for spacer in SpacersLocation['-']:
#     #     countMinus += 1
#     #     for orf in ORFs['-']:
#     #         if (spacer[0] > orf[0]) and (spacer[1] < orf[1]):
#     #             if (Line[spacer[1]-1:spacer[1]-1 + 3] == "CTT") :
#     #                 countInMinus += 1
#     #                 Matches.append([spacer[0],spacer[1], Line[spacer[1] -1 :spacer[1] - 1 + 3], "-"])
#     #                 break
#     #         if ((spacer[0] > orf[0]) and (spacer[1] > orf[1]) and (spacer[0] < orf[1])) or (
#     #                 (spacer[0] < orf[0]) and (spacer[1] < orf[1]) and (spacer[1] > orf[0])):
#     #             if (Line[spacer[1]-1:spacer[1]-1 + 3] == "CTT") :
#     #                 OverlapMinus += 1
#     #                 Matches.append([spacer[0], spacer[1], Line[spacer[1] -1:spacer[1] -1 + 3], "overlapminus"])
#     #                 break
#     #
#     #     for i in range(len(ORFs['-']) - 1):
#     #         if (spacer[0] > ORFs['-'][i][1]) and (spacer[1] < ORFs['-'][i + 1][0]):
#     #             if (Line[spacer[1]-1:spacer[1]-1 + 3] == "CTT") :
#     #                 OutORFminus += 1
#     #                 Matches.append([spacer[0], spacer[1], Line[spacer[1] - 1:spacer[1] - 1 + 3], "outminus"])
#     #                 break
#
#             # else:
#             #     OutORFminus += 1
#             #     Matches.append([spacer[0], spacer[1], Line[spacer[0] - 1:spacer[0] - 1 + 3], "outminus"])
# # with open (OutORFplusFileName, "w") as O:
# #     O.write(OutORFplus)
# # with open (OutORFminusFileName, "w") as O:
# #     O.write(OutORFminus)
# print(countPlus, ' ', countMinus)
# print(countInPlus, ' ' , countInMinus)
# print(OverlapPlus, ' ', OverlapMinus)
# print(OutORFplus,' ', OutORFminus)
# print(InORF, ' ', OutORF)
# with open(MatchesFileName, "w") as MatchesFile:
#      for elem in Matches:
#          MatchesFile.write(str(elem[0]) + ' ' + str(elem[1]) + ' ' + str(elem[2]) + ' ' +str(elem[3]) +' '
#                            +str(elem[4])+ ' ' + str(elem[5])+' ' + str(elem[6]) +"\n")

with open(AAGPositionsFileName, "w") as AAGPositionsFile:
    for elem in AAG_Positions:
        AAGPositionsFile.write(
            str(elem[0]) + ' ' + str(elem[1]) +' ' + str(elem[2]) + "\n")