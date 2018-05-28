import re

EcoliFastaJoinedFileName = "/home/utkinai2/Project1/Ecoli_chromosome_joined.fna"
AAGPositionsFileName = "/home/utkinai2/Project1/AAG_positions_orf.txt"
Naive1MatchesFileName = "/home/utkinai2/Project1/PAMSpacer_matches_AS27first_naive.txt"
Naive2MatchesFileName = "/home/utkinai2/Project1/PAMSpacer_matches_AS27second_naive.txt"
PrimedMatchesFileName = "/home/utkinai2/Project1/PAMSpacer_matches_primed201.txt"

NaivePlusCodingPAMSpacer40FileName = "home/utkinai2/Project1/NaivePlusCodingPAMSpacer4_40nt.txt"
NaiveMinusCodingPAMSpacer40FileName = "home/utkinai2/Project1/NaiveMinusCodingPAMSpacer4_40nt.txt"
PlusCodingPAM40FileName = "home/utkinai2/Project1/PlusCodingPAM_40nt.txt"
MinusCodingPAM40FileName = "home/utkinai2/Project1/MinusCodingPAM_40nt.txt"

for oneLine in open(EcoliFastaJoinedFileName, "r"):
    with open(NaivePlusCodingPAMSpacer40FileName, "w") as NaivePlusCodingPAMSpacer40File:
        with open(NaiveMinusCodingPAMSpacer40FileName, "w") as NaiveMinusCodingPAMSpacer40File:
            for line in open(Naive1MatchesFileName):
                LineValues = line[:-1].split(" ")
                if LineValues[5] == "in":
                    if LineValues[3] == "-":
                        NaiveMinusCodingPAMSpacer40File.write(oneLine[int(line[6]) - 40: int(line[6]) + 40] + "\n")
                    if LineValues[3] == "+":
                        NaivePlusCodingPAMSpacer40File.write(oneLine[int(line[6]) - 40: int(line[6]) + 40] + "\n")
            for line in open(Naive2MatchesFileName):
                LineValues = line[:-1].split(" ")
                if LineValues[5] == "in":
                    if LineValues[3] == "-":
                        NaiveMinusCodingPAMSpacer40File.write(oneLine[int(line[6]) - 40: int(line[6]) + 40] + "\n")
                    if LineValues[3] == "+":
                        NaivePlusCodingPAMSpacer40File.write(oneLine[int(line[6]) - 40: int(line[6]) + 40] + "\n")
