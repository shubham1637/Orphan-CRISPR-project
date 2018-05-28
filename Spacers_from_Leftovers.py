import re
import itertools

SpacersBaseFileName =  "/panfs/pan1.be-md.ncbi.nlm.nih.gov/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/Spacers.fna"
# PartCoveredOrphansFileName = "/home/utkinai2/Project1/Orphan_partcovered.csv"
# PartCoveredSpacersFileName = "/home/utkinai2/Project1/orphans_part_spacers.txt"
# PartCoveredLeftSpacersFileName = "/home/utkinai2/Project1/part_orphans_spacers.txt"
# PartCoveredLeftSpacersFileName1 = "/home/utkinai2/Project1/part_orphans_spacers1.txt"
# PartCoveredLeftSpacersFileName2 = "/home/utkinai2/Project1/part_orphans_spacers2.txt"
# PartCoveredLeftSpacersFileName3 = "/home/utkinai2/Project1/part_orphans_spacers3.txt"
# PartCoveredLeftSpacersFileName4 = "/home/utkinai2/Project1/part_orphans_spacers4.txt"
# PartCoveredLeftSpacersFileName5 = "/home/utkinai2/Project1/part_orphans_spacers5.txt"
# PartCoveredLeftSpacersFileName6 = "/home/utkinai2/Project1/part_orphans_spacers6.txt"
# PartCoveredLeftSpacersFileName7 = "/home/utkinai2/Project1/part_orphans_spacers7.txt"
# PartCoveredLeftSpacersFileName8 = "/home/utkinai2/Project1/part_orphans_spacers8.txt"
# FullCoveredOrphansFileName = "/home/utkinai2/Project1/Orphan_goodFullcovered.csv"
# FullCoveredSpacersFileName = "/home/utkinai2/Project1/orphans_full_spacers.txt"
# FullCoveredLeftSpacersFileName = "/home/utkinai2/Project1/full_orphans_spacers.txt"
# FullCoveredLeftSpacersFileName1 = "/home/utkinai2/Project1/full_orphans_spacers1.txt"
# FullCoveredLeftSpacersFileName2 = "/home/utkinai2/Project1/full_orphans_spacers2.txt"
# FullCoveredLeftSpacersFileName3 = "/home/utkinai2/Project1/full_orphans_spacers3.txt"
# FullCoveredLeftSpacersFileName4 = "/home/utkinai2/Project1/full_orphans_spacers4.txt"
# FullCoveredLeftSpacersFileName5 = "/home/utkinai2/Project1/full_orphans_spacers5.txt"
# FullCoveredLeftSpacersFileName6 = "/home/utkinai2/Project1/full_orphans_spacers6.txt"
# FullCoveredLeftSpacersFileName7 = "/home/utkinai2/Project1/full_orphans_spacers7.txt"
# FullCoveredLeftSpacersFileName8 = "/home/utkinai2/Project1/full_orphans_spacers8.txt"

OrphansWithProteinHitsFileName = "/home/utkinai2/Project1/Orphans_withProtHits.txt"
SpacersFromProteinHitsFileName = "/home/utkinai2/Project1/Orphans_withProtHits_spacers.txt"
OrphansWithSameSpacersFileName = "/home/utkinai2/Project1/Orphans_with_sameSpacers.txt"
SpacersFromBadOrphansFileName = "/home/utkinai2/Project1/Orphans_with_sameSpacers_spacers.txt"

Orphans99FileName = "/home/utkinai2/Project1/99newOrphans.txt"
SpacersFrom99OrphansFileName = "/home/utkinai2/Project1/99newOrphans_spacers.txt"

# BadLeftSpacersFileName = "/home/utkinai2/Project1/bad_orphans_spacers.txt"
# BadLeftSpacersFileName1 = "/home/utkinai2/Project1/bad_orphans_spacers1.txt"
# BadLeftSpacersFileName2 = "/home/utkinai2/Project1/bad_orphans_spacers2.txt"
# BadLeftSpacersFileName3 = "/home/utkinai2/Project1/bad_orphans_spacers3.txt"
# BadLeftSpacersFileName4 = "/home/utkinai2/Project1/bad_orphans_spacers4.txt"
# BadLeftSpacersFileName5 = "/home/utkinai2/Project1/bad_orphans_spacers5.txt"
# BadLeftSpacersFileName6 = "/home/utkinai2/Project1/bad_orphans_spacers6.txt"
# BadLeftSpacersFileName7 = "/home/utkinai2/Project1/bad_orphans_spacers7.txt"
# BadLeftSpacersFileName8 = "/home/utkinai2/Project1/bad_orphans_spacers8.txt"

count = 0
Const = 2000
SpacerID = {}
#creation of a dict, where the ID is a key and CrisprStart-CrisprEnd is a value (boundaries to search in a Spacers.fna file)
for line in open(Orphans99FileName, "r"):
    count += 1
    if count < 2:
        continue  # header is skipped now
    LineValues = line[:-1].split(" ")
    IDMatch = re.search(r'([A-Z]{2,}([0-9]{1,}))', LineValues[2])
    if IDMatch:
        # print(str(IDMatch.group(0)))
       #SpacerID[str(IDMatch.group(0))] = [LineValues[6],LineValues[7]]
        SpacerID[str(IDMatch.group(0))] = [LineValues[3], LineValues[4]]

#search
with open(SpacersFrom99OrphansFileName, "w") as FullCoveredSpacersFile:
    with open(SpacersBaseFileName, "r") as f:
        for line1, line2 in itertools.zip_longest(*[f]* 2):
            LineValues = line1[:-1].split("_")
            IDMatch = re.search(r'([A-Z]{2,}([0-9]{1,}))', LineValues[0])
            if IDMatch:
                if str(IDMatch.group(0)) in SpacerID.keys():
                    FullCoveredSpacersFile.write(line1 + line2)

        # if LineValues[0] == ">"
        # if LineValues[2].strip("\"") in ContigIDSizes:
        #     if (int(LineValues[4]) < Offset) or (int(LineValues[5]) > ContigIDSizes[LineValues[2].strip("\"")] - Offset):
        #         continue
        #     LeftCRISPRsWithinContigFile.write(line)

# with open(Spacers_FromProteinHitsFileName, "w") as FullCoveredLeftSpacersFile:
#     for line in open(SpacersFromProteinHitsFileName, "r"):
#         LineValues = line[:-1].split(" ")
#         if LineValues[0] == " ":
#             continue
#         FullCoveredLeftSpacersFile.write(line)
count1 = 0
# with open(BadLeftSpacersFileName1, "w") as PartCoveredLeftSpacersFile1:
#     with open(BadLeftSpacersFileName2, "w") as PartCoveredLeftSpacersFile2:
#         with open(BadLeftSpacersFileName3, "w") as PartCoveredLeftSpacersFile3:
#             with open(BadLeftSpacersFileName4, "w") as PartCoveredLeftSpacersFile4:
#                 with open(BadLeftSpacersFileName5, "w") as PartCoveredLeftSpacersFile5:
#                     with open(BadLeftSpacersFileName6, "w") as PartCoveredLeftSpacersFile6:
#                         with open(BadLeftSpacersFileName7, "w") as PartCoveredLeftSpacersFile7:
#                             with open(BadLeftSpacersFileName8, "w") as PartCoveredLeftSpacersFile8:
#                                 for line in open(SpacersFromBadOrphansFileName, "r"):
#                                     count1 += 1
#                                     if (count1 > 0) and (count1 <= 32000):
#                                         PartCoveredLeftSpacersFile1.write(line)
#                                     if (count1 > 32000) and (count1 <= 64000):
#                                         PartCoveredLeftSpacersFile2.write(line)
#                                     if (count1 > 64000) and (count1 <= 96000):
#                                         PartCoveredLeftSpacersFile3.write(line)
#                                     if (count1 > 96000) and (count1 <= 128000):
#                                         PartCoveredLeftSpacersFile4.write(line)
#                                     if (count1 > 128000) and (count1 <= 160000):
#                                         PartCoveredLeftSpacersFile5.write(line)
#                                     if (count1 > 160000) and (count1 <= 192000):
#                                         PartCoveredLeftSpacersFile6.write(line)
#                                     if (count1 > 192000) and (count1 <= 225000):
#                                         PartCoveredLeftSpacersFile7.write(line)
#                                     if (count1 > 225000) and (count1 <= 255000):
#                                         PartCoveredLeftSpacersFile8.write(line)








