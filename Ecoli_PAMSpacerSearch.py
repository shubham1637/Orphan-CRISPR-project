import re

UniquePAMSpacersFileName = "//frosty/utkinai2/Project1/PAMSpacer_matches_primed201.txt"
EcoliORFsFileName = "C:/Users/utkinai2/Desktop/Ira/Ecoli_chromosome.pty"

PositionPlusAAG = []
PositionMinusAAG = []
PositionPlusCTT = []
PositionMinusCTT=[]
count = 0

ORFs = {'+': [], '-': []}
ORFs2 = []
for Line in open(EcoliORFsFileName, "r"):
    if Line[0] == "#":
        continue
    LineValues = Line[:-1].split("\t")  # the last one is perenos stroki, thats'why -1
    RangeMatch = re.search(r'([0-9]{1,}..[0-9]{1,})', LineValues[1])
    Range = RangeMatch.group(0).split("..")

    ORFs2.append([int(Range[0]), int(Range[1])])

    if LineValues[2] == "+":
        ORFs['+'].append([int(Range[0]), int(Range[1])])
    if LineValues[2] == "-":
        ORFs['-'].append([int(Range[0]), int(Range[1])])

for line in open(UniquePAMSpacersFileName, "r"):
    count += 1
    if count < 2:
        continue
    LineValues = line[:-1].split(" ")
    # for orf in ORFs2:
    #      if (int(LineValues[0]) > orf[0]) and (int(LineValues[1]) < orf[1]):
    #               PositionPlusAAG.append((int(LineValues[3]) - orf[0] ) % 3 +1)
    if LineValues[3] == "+" and LineValues[4] == "AAG":
        for orf in ORFs['+']:
            if (int(LineValues[0]) > orf[0]) and (int(LineValues[1]) < orf[1]):
                PositionPlusAAG.append((int(LineValues[0]) - orf[0]) % 3 + 1)
    if LineValues[3] == "-" and LineValues[4] == "AAG":
            for orf in ORFs['-']:
               if (int(LineValues[0]) > orf[0]) and (int(LineValues[1]) < orf[1]):
                  PositionMinusAAG.append((orf[1] - int(LineValues[1]) ) % 3 + 1)

    if LineValues[3] == "+" and LineValues[4] == "CTT":
            for orf in ORFs['+']:
               if (int(LineValues[0]) > orf[0]) and (int(LineValues[1]) < orf[1]):
                  PositionPlusCTT.append((int(LineValues[0]) - orf[0] ) % 3 + 1)
    if LineValues[3] == "-" and LineValues[4] == "CTT":
            for orf in ORFs['-']:
               if (int(LineValues[0]) > orf[0]) and (int(LineValues[1]) < orf[1]):
                  PositionMinusCTT.append((int(LineValues[0]) - orf[0] ) % 3 +1)

# for match in re.finditer('AAG', 'AAGGRTEGDAAGAAUJSJSAAGHJEEJ'):
#       PositionPlusAAG.append([match.start(),(int(match.start()) + 3) % 3])
print(PositionPlusAAG)

for count, elem in sorted(((PositionPlusAAG.count(e), e) for e in set(PositionPlusAAG)), reverse=True):
    print ('PositionPlusAAG %s (%d)' % (elem, count))
for count, elem in sorted(((PositionMinusAAG.count(e), e) for e in set(PositionMinusAAG)), reverse=True):
    print ('PositionMinusAAG %s (%d)' % (elem, count))
for count, elem in sorted(((PositionPlusCTT.count(e), e) for e in set(PositionPlusCTT)), reverse=True):
    print ('PositionPlusCTT %s (%d)' % (elem, count))
for count, elem in sorted(((PositionMinusCTT.count(e), e) for e in set(PositionMinusCTT)), reverse=True):
    print ('PositionMinusCTT %s (%d)' % (elem, count))