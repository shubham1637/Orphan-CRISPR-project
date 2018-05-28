import re
import subprocess

EcoliORFsFileName = "/home/utkinai2/Project1/Ecoli_chromosome.pty"
EcoliFastaFileName = "/home/utkinai2/Project1/Ecoli_chromosome.fna"
EcoliFastaJoinedFileName = "/home/utkinai2/Project1/Ecoli_chromosome_joined.fna"
EcoliPlusStrandORFsFileName = "/home/utkinai2/Project1/Ecoli_plusstrand.txt"
EcoliMinusStrandORFsFileName = "/home/utkinai2/Project1/Ecoli_minusstrand.txt"

EcoliMinusFileName = "/home/utkinai2/Project1/Ecoli_minus.txt"
EcoliNonCodingFileName  = "/home/utkinai2/Project1/Ecoli_non_coding.txt"

count = 0
PlusCount = 0
Positions_G_PlusStrand = []
Positions_afterG_PlusStrand = []

with open(EcoliFastaJoinedFileName, "w") as EcoliFastaJoinedFile:
    with open(EcoliFastaFileName, "r") as f:
        EcoliFastaJoinedFile.write("".join(line.strip() for line in f))

Ranges = []
for Line in open(EcoliORFsFileName, "r"):
    count += 1
    if count < 2:
        continue  # header is skipped now
    LineValues = Line[:-1].split("\t") # the last one is perenos stroki, thats'why -1
    RangeMatch = re.search(r'([0-9]{1,}..[0-9]{1,})', LineValues[1])
    Range = RangeMatch.group(0).split("..")
    if LineValues[2] == "+":
        Ranges.append([Range[0], Range[1], "+"])
    if LineValues[2] == "-":
        Ranges.append([Range[0], Range[1], "-"])

for Line in open(EcoliFastaJoinedFileName, "r"):
    for elem in Ranges:
        if elem[2] == "+":
            if re.search(r'AAG', Line[elem[0]:elem[1]]):
                if int(re.search(r'AAG', elem).start()) + 1 + 33 <= elem[1]:
                            PlusCount += 1
                            Positions_G_PlusStrand.append((int(re.search(r'AAG', elem).start()) + 2) % 3)
                            Positions_afterG_PlusStrand.append((int(re.search(r'AAG', elem).start()) + 3) % 3)
            if int(re.search(r'AAG', elem).start()) + 2 + 33 <= len(elem):
                Positions_afterG_PlusStrand.append((int(re.search(r'AAG', elem).start()) + 3) % 3)
            else:
                CountPlus += 1
print('OUT OF PROTEIN (plus) ', CountPlus, ' times')