import re

SourceFileName = "/home/utkinai2/Project1/merged_all_with_reverserepeats.csv"
#RepeatsFileName = "/home/utkinai2/Project1/identified_filter.txt"
RepeatsAsFastaFileName = "/home/utkinai2/Project1/identified_orphans_repeats_withreverse.txt"

count = 0
IDandRange = {}
# for line in open(SourceFileName, "r"):
#     count += 1
#     if count < 2:
#         continue  # header is skipped now
#     LineValues = line[:-1].split(",")
    #RangeMatch = re.search(r'([0-9]{1,}\'\,\s\'[0-9]{1,})', LineValues[1])
    #if RangeMatch:
#     IDandRange[LineValues[0]] = [int(RangeMatch.group(0).split("', '")[0]), int(RangeMatch.group(0).split("', '")[1])]
# print(IDandRange)


with open(RepeatsAsFastaFileName, "w") as RepeatsAsFastaFile:
    for line in open(SourceFileName, "r"):
        # count += 1
        # if count < 2:
        #     continue  # header is skipped now
        LineValues = line[:-1].split(",")
        #ID = re.search( r'([A-Z]{1,}[0-9]{1,})', LineValues[1]).group(0)
        RepeatsAsFastaFile.write(">" + str(LineValues[0])+ '\n' + LineValues[6]+'\n')