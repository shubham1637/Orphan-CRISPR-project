import re

#make Reverse complement of EcoliMinusStrandORFFileName and
ReverseMinusStrandFileName = "/home/utkinai2/Project1/Ecoli_minusstrand_reversecompl.txt"
ReverseMinusFormattedFileName =  "/home/utkinai2/Project1/Ecoli_minusstrand_reverseform.txt"
EcoliPlusStrandORFsFileName = "/home/utkinai2/Project1/Ecoli_plusstrand.txt"
EcoliPlusStrandFormattedFileName = "/home/utkinai2/Project1/Ecoli_plusstrand_form.txt"
EcoliNonCodingFileName = "/home/utkinai2/Project1/Ecoli_noncoding.txt"
EcoliNonCodingFormattedFileName = "/home/utkinai2/Project1/Ecoli_noncoding_form.txt"
EcoliFastaJoinedFileName = "/home/utkinai2/Project1/Ecoli_chromosome_joined.fna"
# EcoliAllFileName = "/home/utkinai2/Project1/Ecoli_all.txt"
# EcoliAllFormattedFileName = "/home/utkinai2/Project1/Ecoli_all_form.txt"

with open(ReverseMinusFormattedFileName, "w") as ReverseMinusFormattedFile:
    with open(ReverseMinusStrandFileName, "r") as f:
        ReverseMinusFormattedFile.write("".join(line.strip() for line in f))

with open(EcoliPlusStrandFormattedFileName, "w") as EcoliPlusStrandFormattedFile:
    with open(EcoliPlusStrandORFsFileName, "r") as f:
        EcoliPlusStrandFormattedFile.write("".join(line.strip() for line in f))

with open(EcoliNonCodingFormattedFileName, "w") as EcoliNonCodingFormattedFile:
    with open(EcoliNonCodingFileName, "r") as f:
        EcoliNonCodingFormattedFile.write("".join(line.strip() for line in f))

# with open(EcoliAllFormattedFileName, "w") as EcoliAllFormattedFile:
#     with open(EcoliAllFileName, "r") as f:
#         EcoliAllFormattedFile.write("".join(line.strip() for line in f))

CountMinusAAG = 0
CountPlusAAG = 0
CountMinusCTT = 0
CountPlusCTT = 0
CountPlusTTC = 0
CountMinusTTC = 0
CountPlusGAA = 0
CountMinusGAA = 0
CountPlusCAG = 0
CountMinusCAG = 0

CountNonCodingAAG = 0
CountNonCodingCTT = 0
CountNonCodingGAA = 0
CountNonCodingTTC = 0

OutPlusAAG = 0
OutMinusAAG = 0
OutPlusCTT = 0
OutMinusCTT =0
OutPlusTTC = 0
OutMinusTTC =0
OutPlusGAA =0
OutMinusGAA = 0
OutMinusCAG = 0
OutPlusCAG = 0


Positions_G_PlusStrandAAG = []
Positions_afterG_PlusStrandAAG = []
Positions_G_PlusStrandCAG = []
Positions_afterG_PlusStrandCAG = []
Positions_G_PlusStrandCTT = []
Positions_afterG_PlusStrandCTT = []
Positions_G_PlusStrandTTC = []
Positions_afterG_PlusStrandTTC = []
Positions_G_PlusStrandGAA = []
Positions_afterG_PlusStrandGAA = []

for Line in open(EcoliPlusStrandFormattedFileName, "r"):
    LineValues = Line[:-1].split(">CP009273.1 Escherichia coli BW25113, complete genome")
    for elem in LineValues:
        if re.search(r'AAG',elem):
            CountPlusAAG += elem.count("AAG")
            for match in re.finditer(r'AAG', elem):
                if int(match.start())+ 1 + 33 <= len(elem):
                    Positions_G_PlusStrandAAG.append((int(match.start()) + 2) % 3 + 1)
                if int(match.start()) + 2 + 33 <= len(elem):
                    Positions_afterG_PlusStrandAAG.append((int(match.start()) + 3) % 3 +1)
                else:
                    OutPlusAAG +=1
        if re.search(r'CAG', elem):
            CountPlusCAG += elem.count("CAG")
            for match in re.finditer(r'CAG', elem):
                if int(match.start()) + 1 + 33 <= len(elem):
                    Positions_G_PlusStrandCAG.append((int(match.start()) + 2) % 3 + 1 )
                if int(match.start()) + 2 + 33 <= len(elem):
                        Positions_afterG_PlusStrandCAG.append((int(match.start()) + 3) % 3 + 1)
                else:
                    OutPlusCAG += 1
        # if re.search(r'CTT', elem):
        #     CountPlusCTT += elem.count("CTT")
        #     for match in re.finditer(r'CTT', elem):
        #         if int(match.start()) + 1 + 33 <= len(elem):
        #             Positions_G_PlusStrandCTT.append((int(match.start()) + 2) % 3 + 1)
        #         if int(match.start()) + 2 + 33 <= len(elem):
        #             Positions_afterG_PlusStrandCTT.append((int(match.start()) + 3) % 3 + 1)
        #         else:
        #             OutPlusCTT += 1
        # if re.search(r'TTC', gene):
        #     CountPlusTTC += gene.count("TTC")
        #     for match in re.finditer(r'TTC', gene):
        #         if int(match.start()) + 1 + 33 <= len(gene):
        #             Positions_G_PlusStrandTTC.append((int(match.start()) + 2) % 3)
        #         if int(match.start()) + 2 + 33 <= len(gene):
        #             Positions_afterG_PlusStrandTTC.append((int(match.start()) + 3) % 3)
        #         else:
        #             OutPlusTTC += 1
        # if re.search(r'GAA', gene):
        #     CountPlusGAA += gene.count("GAA")
        #     for match in re.finditer(r'GAA', gene):
        #         if int(match.start()) + 1 + 33 <= len(gene):
        #             Positions_G_PlusStrandGAA.append((int(match.start()) + 2) % 3)
        #         if int(match.start()) + 2 + 33 <= len(gene):
        #             Positions_afterG_PlusStrandGAA.append((int(match.start()) + 3) % 3)
        #         else:
        #             OutPlusGAA += 1
print('OUT OF PROTEIN AAG (plus) ', OutPlusAAG, ' times')
print('OUT OF PROTEIN CAG (plus) ', OutPlusCAG, ' times')

CountAAG = 0
for line in open(EcoliFastaJoinedFileName, "r"):
    for match in re.finditer(r'AAG', line):
        CountAAG += 1

print(CountAAG)

Positions_G_MinusStrandAAG = []
Positions_afterG_MinusStrandAAG = []
Positions_G_MinusStrandCAG = []
Positions_afterG_MinusStrandCAG = []
Positions_G_MinusStrandCTT = []
Positions_afterG_MinusStrandCTT = []
Positions_G_MinusStrandTTC = []
Positions_afterG_MinusStrandTTC = []
Positions_G_MinusStrandGAA = []
Positions_afterG_MinusStrandGAA = []

for Line in open(ReverseMinusFormattedFileName, "r"):
    LineValues = Line[:-1].split(">CP009273.1 Escherichia coli BW25113, complete genome reverse complement")
    for elem in LineValues:
        if re.search(r'AAG',elem):
            CountMinusAAG += elem.count("AAG")
            for match in re.finditer(r'AAG', elem):
                if int(match.start()) - 33 >= 0:
                    Positions_G_MinusStrandAAG.append((int(match.start()) + 2) % 3 + 1)
                # if int(match.start()) + 1 - 33 >= 0:
                #     Positions_afterG_MinusStrandAAG.append((int(match.start()) + 3) % 3 + 1)
                else:
                     OutMinusAAG += 1
        if re.search(r'CAG', elem):
            CountMinusCAG += elem.count("CAG")
            for match in re.finditer(r'CAG', elem):
                if int(match.start()) - 33 >= 0:
                    Positions_G_MinusStrandCAG.append((int(match.start()) + 2) % 3 +1 )
                # if int(match.start()) + 1 - 33 >= 0:
                #     Positions_afterG_MinusStrandCAG.append((int(match.start()) + 3) % 3 + 1 )
                else:
                    OutMinusCAG += 1
        # if re.search(r'CTT',elem):
        #     CountMinusCTT += elem.count("CTT")
        #     for match in re.finditer(r'CTT', elem):
        #         if int(match.start()) - 33 >= 0:
        #             Positions_G_MinusStrandCTT.append((int(match.start())+ 2) % 3 + 1)
        #         if int(match.start()) + 1 - 33 >= 0:
        #             Positions_afterG_MinusStrandCTT.append((int(match.start()) + 3) % 3 +1 )
        #         else:
        #              OutMinusCTT += 1
        # if re.search(r'TTC', gene):
        #     CountMinusTTC += gene.count("TTC")
        #     for match in re.finditer(r'TTC', gene):
        #         if int(match.start()) + 1 + 33 <= len(gene):
        #             Positions_G_MinusStrandTTC.append((int(match.start()) + 2) % 3)
        #         if int(match.start()) + 2 + 33 <= len(gene):
        #             Positions_afterG_MinusStrandTTC.append((int(match.start()) + 3) % 3)
        #         else:
        #             OutMinusTTC += 1
        # if re.search(r'GAA', gene):
        #     CountMinusGAA += gene.count("GAA")
        #     for match in re.finditer(r'GAA', gene):
        #         if int(match.start()) + 1 + 33 <= len(gene):
        #             Positions_G_MinusStrandGAA.append((int(match.start()) + 2) % 3)
        #         if int(match.start()) + 2 + 33 <= len(gene):
        #             Positions_afterG_MinusStrandGAA.append((int(match.start()) + 3) % 3)
        #         else:
        #             OutMinusGAA += 1
print('OUT OF PROTEIN AAG (minus) ', OutMinusAAG, ' times')
print('OUT OF PROTEIN CAG (minus) ', OutMinusCAG, ' times')

for Line in open(EcoliNonCodingFormattedFileName, "r"):
    LineValues = Line[:-1].split(">CP009273.1 Escherichia coli BW25113, complete genome")
    for elem in LineValues:
        if re.search(r'AAG',elem):
            CountNonCodingAAG += elem.count("AAG")
        if re.search(r'CTT',elem):
            CountNonCodingCTT += elem.count("CTT")
        if re.search(r'TTC', elem):
            CountNonCodingTTC += elem.count("AAG")
        if re.search(r'GAA', elem):
            CountNonCodingGAA += elem.count("CTT")
print('PAM AAG in NonCoding region ', CountNonCodingAAG, ' times')
print('PAM CTT in NonCoding region ', CountNonCodingCTT, ' times')
print('PAM GAA in NonCoding region ', CountNonCodingGAA, ' times')
print('PAM TTC in NonCoding region ', CountNonCodingTTC, ' times')


# for Line in open(EcoliAllFormattedFileName, "r"):
#     LineValues = Line[:-1].split(">CP009273.1 Escherichia coli BW25113, complete genome")
#     for gene in LineValues:
#         if re.search(r'AAG',gene):
#             CountAll += 1
# print('PAM in all genome ', CountAll, ' times')

for count, elem in sorted(((Positions_G_PlusStrandAAG.count(e), e) for e in set(Positions_G_PlusStrandAAG)), reverse=True):
    print ('Positions_G_PlusStrandAAG %s (%d)' % (elem, count))
# for count, elem in sorted(((Positions_afterG_PlusStrandAAG.count(e), e) for e in set(Positions_afterG_PlusStrandAAG)), reverse=True):
#     print ('Positions_afterG_PlusStrandAAG %s (%d)' % (elem, count))
for count, elem in sorted(((Positions_G_MinusStrandAAG.count(e), e) for e in set(Positions_G_MinusStrandAAG)), reverse=True):
    print ('Positions_G_MinusStrandAAG %s (%d)' % (elem, count))
# for count, elem in sorted(((Positions_afterG_MinusStrandAAG.count(e), e) for e in set(Positions_afterG_MinusStrandAAG)), reverse=True):
#     print ('Positions_afterG_MinusStrandAAG %s (%d)' % (elem, count))
for count, gene in sorted(((Positions_G_PlusStrandCAG.count(e), e) for e in set(Positions_G_PlusStrandCAG)), reverse=True):
    print ('Positions_G_PlusStrandCAG %s (%d)' % (gene, count))
# for count, gene in sorted(((Positions_afterG_PlusStrandCAG.count(e), e) for e in set(Positions_afterG_PlusStrandCAG)), reverse=True):
#     print ('Positions_afterG_PlusStrandCAG %s (%d)' % (gene, count))
for count, gene in sorted(((Positions_G_MinusStrandCAG.count(e), e) for e in set(Positions_G_MinusStrandCAG)), reverse=True):
    print ('Positions_G_MinusStrandCAG %s (%d)' % (gene, count))
# for count, gene in sorted(((Positions_afterG_MinusStrandCAG.count(e), e) for e in set(Positions_afterG_MinusStrandCAG)), reverse=True):
#     print ('Positions_afterG_MinusStrandCAG %s (%d)' % (gene, count))
for count, elem in sorted(((Positions_G_PlusStrandCTT.count(e), e) for e in set(Positions_G_PlusStrandCTT)), reverse=True):
    print ('Positions_G_PlusStrandCTT %s (%d)' % (elem, count))
for count, elem in sorted(((Positions_afterG_PlusStrandCTT.count(e), e) for e in set(Positions_afterG_PlusStrandCTT)), reverse=True):
    print ('Positions_afterG_PlusStrandCTT %s (%d)' % (elem, count))
for count, elem in sorted(((Positions_G_MinusStrandCTT.count(e), e) for e in set(Positions_G_MinusStrandCTT)), reverse=True):
    print ('Positions_G_MinusStrandCTT %s (%d)' % (elem, count))
# for count, elem in sorted(((Positions_afterG_MinusStrandCTT.count(e), e) for e in set(Positions_afterG_MinusStrandCTT)), reverse=True):
#     print ('Positions_afterG_MinusStrandCTT %s (%d)' % (elem, count))
# for count, gene in sorted(((Positions_G_PlusStrandGAA.count(e), e) for e in set(Positions_G_PlusStrandGAA)), reverse=True):
#     print ('Positions_G_PlusStrandGAA %s (%d)' % (gene, count))
# for count, gene in sorted(((Positions_afterG_PlusStrandGAA.count(e), e) for e in set(Positions_afterG_PlusStrandGAA)), reverse=True):
#     print ('Positions_afterG_PlusStrandGAA %s (%d)' % (gene, count))
# for count, gene in sorted(((Positions_G_MinusStrandGAA.count(e), e) for e in set(Positions_G_MinusStrandGAA)), reverse=True):
#     print ('Positions_G_MinusStrandGAA %s (%d)' % (gene, count))
# for count, gene in sorted(((Positions_afterG_MinusStrandGAA.count(e), e) for e in set(Positions_afterG_MinusStrandGAA)), reverse=True):
#     print ('Positions_afterG_MinusStrandGAA %s (%d)' % (gene, count))
# for count, gene in sorted(((Positions_G_PlusStrandTTC.count(e), e) for e in set(Positions_G_PlusStrandTTC)), reverse=True):
#     print ('Positions_G_PlusStrandTTC %s (%d)' % (gene, count))
# for count, gene in sorted(((Positions_afterG_PlusStrandTTC.count(e), e) for e in set(Positions_afterG_PlusStrandTTC)), reverse=True):
#     print ('Positions_afterG_PlusStrandTTC %s (%d)' % (gene, count))
# for count, gene in sorted(((Positions_G_MinusStrandTTC.count(e), e) for e in set(Positions_G_MinusStrandTTC)), reverse=True):
#     print ('Positions_G_MinusStrandTTC %s (%d)' % (gene, count))
# for count, gene in sorted(((Positions_afterG_MinusStrandTTC.count(e), e) for e in set(Positions_afterG_MinusStrandTTC)), reverse=True):
#     print ('Positions_afterG_MinusStrandTTC %s (%d)' % (gene, count))
#
