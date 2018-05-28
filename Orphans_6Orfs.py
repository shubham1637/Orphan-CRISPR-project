import subprocess

OrphansWithORFFileName = "/home/utkinai2/Project1/identified_translated_in6ORFs.txt"
OrphansFileName = "/home/utkinai2/Project1/identified.csv"
ORFResultsFileName = "Identified_6orfs.txt"
#FinalFileName = "/home/utkinai2/Project1/orphans_translated_in6ORFs_.txt"

# if orf doesn't work - export PATH=$PATH:~wolf/bin/ and then do: export PATH=$PATH:~wolf/bin/add/

count = 0
count1 = 0
with open(OrphansWithORFFileName, "w") as OrphansWithORFFile:
    for Line in open(OrphansFileName, "r"):
        count += 1
        if count < 2:
            continue  # header is skipped now
        LineValues = Line[:-1].split(",")  # the last one is perenos stroki, thats'why -1
        subprocess.call("blastdbcmd -db /panfs/pan1/prokdata/db/all1603.nt" + " -entry " + LineValues[2] +
                        " -range " + str(LineValues[3]) + "-" + str(LineValues[4]) + " | orf -f=0 -m=1 -n=1 -l=40 > " + ORFResultsFileName, shell=True)
        for line in open(ORFResultsFileName, "r"):
            if line[0] == ">":
                count1 += 1
                LineVal = line[:-1].split(" ")
                OrphansWithORFFile.write('>' + LineVal[1] + '_' + str(LineValues[3]) + '_' + str(LineValues[4]) + '_' + str(count1) + '\t' +
                                     LineVal[-6] + '\t' + LineVal[-5] + '\t' + LineVal[-4] + '\t' + LineVal[-3] + '\t' + LineVal[-2] + '\t' + LineVal[-1]+ '\n')
            else:
                OrphansWithORFFile.write(line)

# with open(FinalFileName, "w") as FinalFile:
#     for line in open(OrphansWithORFFileName, "r"):
#         if line[0] == ">":
#             count1 += 1
#             LineVal = line[:-1].split("\t")
#             FinalFile.write(LineVal[0] + '_' + str(count1) + '\t' + LineVal[1] + '\t' + LineVal[2] + '\t' + LineVal[3] + '\t' + LineVal[4] + '\t' + LineVal[5] + '\t' + LineVal[6] + '\n')
#         else:
#             FinalFile.write(line)
