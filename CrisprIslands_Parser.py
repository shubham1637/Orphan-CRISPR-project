import re

LeftoversFilterFileName = "/home/utkinai2/Project1/identified_filter.txt"
IslandsAnnotationFileName = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/prokdata/CRISPRicity/SpacerAnalysis/AllIteration/TmpFolder/CRISPR_islands_new/Islands.ann_CRISPR"
LeftoversNearCasFileName = "/home/utkinai2/Project1/identified_filter_near_cas1.txt"

LeftoversIDandRange = {}
count = 0
Const = 2000

for line in open(LeftoversFilterFileName, "r"):
    # count += 1
    # if count < 2:
    #     continue  # header is skipped now
    LineValues = line[:-1].split(",")
    IDMatch = re.search(r'([A-Z]{2,}([0-9]{1,}))', LineValues[2])
    if IDMatch:
        # print(str(IDMatch.group(0)))
        LeftoversIDandRange[str(IDMatch.group(0))] = [LineValues[3],LineValues[4]]




with open(LeftoversNearCasFileName, "w") as LeftoversNearCasFile:
    for line in open(IslandsAnnotationFileName, "r"):
        LineValues = line[:-1].split("\t")  # в строках-заголовках [1] - это пробел!!! а в строках обычных [7] - это пробел!!! (wtf, but keep in mind)
        if LineValues[0] == "===":
            Numbers = [str(i) for i in str(LineValues[1]).split(",")] # list with numbers of rows, where are cas genes
            #print(Numbers)
        else:
            if (Numbers[0].isdigit() == True) and (LineValues[4] in LeftoversIDandRange.keys()):
                   # если в списке не пробелы, а действительно номера строк, и рассматривая строка в этом списке, и ID есть в словаре Лефтоверов
                if LineValues[7] == "+":
                # RangeMatch = re.search(r'([0-9]{1,})..([0-9]{1,})', LineValues[1])
                # if RangeMatch:
                #     range = str(RangeMatch.group(0))
                #     range = [int(i) for i in range.split("..")]
                # if ((range[0] - Const) < int(LeftoversIDandRange[str(LineValues[4])][0])) or \
                #         ((range[1] + Const) > int(LeftoversIDandRange[str(LineValues[4])][1])):
                # если начало криспр больше чем на 10кб дальше, чем начало рамки из файла CRISPRislands или конец аналогично

                    if  (int(LeftoversIDandRange[LineValues[4]][0]) > int(int(re.search(r'([0-9]{1,})..([0-9]{1,})', LineValues[1]).group(0).split("..")[0]) - 10000) and
                             (int(LeftoversIDandRange[LineValues[4]][0]) < int(int(re.search(r'([0-9]{1,})..([0-9]{1,})', LineValues[1]).group(0).split("..")[1])) + 10000) ) or \
                            ( (int(LeftoversIDandRange[LineValues[4]][1]) > int(re.search(r'([0-9]{1,})..([0-9]{1,})', LineValues[1]).group(0).split("..")[0]) - 10000) and
                                  (int(LeftoversIDandRange[LineValues[4]][1]) < int(re.search(r'([0-9]{1,})..([0-9]{1,})', LineValues[1]).group(0).split("..")[1])+ 10000) ):
                     #if re.search(r'(cas[0-9]{1,}|csx|csm)', LineValues[8]):
                     if re.search(r'((?:^|\W)cas1(?:$|\W))', LineValues[8]):
                        LeftoversNearCasFile.write(str(LineValues[4]) + ' ' + str(LeftoversIDandRange[LineValues[4]]) + ' '+ str(LineValues[1])+' ' + str(LineValues[8]) + '\n')
    #     else:
    #         continue
    # if c == 256:
    #     break






