

PhylumFileName = "/panfs/pan1/prokdata/db_tmp/all1603.tax.tab"

count = 0
Last = ''
for line in open(PhylumFileName, "r"):
    LineValues = line[:-1].split('\t')
    if count > 50:
        break

    if LineValues[-1] != Last:
        print(LineValues[0], ' BEGINNING')
        print(LineValues[-1])
        print(LineValues[-2])
        print(LineValues[-3])
        print(LineValues[-4])
    Last = LineValues[-1]