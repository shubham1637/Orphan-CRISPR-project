import re
import subprocess
import itertools

EcoliORFsFileName = "/home/utkinai2/Project1/Ecoli_chromosome.pty"
EcoliNonCodingFileName = "/home/utkinai2/Project1/Ecoli_noncoding.txt"
EcoliAllFileName = "/home/utkinai2/Project1/Ecoli_all.txt"

count = 0
Ranges = []

for line in open(EcoliORFsFileName, "r"):
    count += 1
    if count < 2:
        continue  # header is skipped now
    LineValues = line[:-1].split("\t")  # the last one is perenos stroki, thats'why -1
    RangeMatch = re.search(r'([0-9]{1,}..[0-9]{1,})', LineValues[1])
    Range = RangeMatch.group(0).split("..")
    Ranges.append(Range[0])
    Ranges.append(Range[1])
    if count > 3:
    # if LineValues[2] == "+":
    #     plusProteins.append([LineValues[5],LineValues[6]])
        subprocess.call("blastdbcmd -db /panfs/pan1/prokdata/db/all1603.nt" + " -entry " + LineValues[4] +
                     " -range " + str(int(Ranges[-3])) + "-" + str(int(Range[-2])) + " >> " + EcoliNonCodingFileName, shell=True)
        subprocess.call("blastdbcmd -db /panfs/pan1/prokdata/db/all1603.nt" + " -entry " + LineValues[4] +
                    " -range " + str(int(Ranges[0])) + "-" + str(int(Range[-1])) + " >> " + EcoliAllFileName, shell=True)
    # if LineValues[2] == "-":
    #     minusProteins.append([LineValues[5], LineValues[6]])
    #     subprocess.call("blastdbcmd -db /panfs/pan1/prokdata/db/all1603.nt" + " -entry " + LineValues[4] +
    #                 " -range " + str(int(Range[0])) + "-" + str(int(Range[1])) + " >> " + EcoliMinusStrandORFsFileName, shell=True)
