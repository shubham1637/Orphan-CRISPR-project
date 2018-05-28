import io
import itertools
import re


PartCoveredSpacersFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/part_orphans_spacers.txt"
GoodFullCoveredSpacersFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/full_orphans_spacers.txt"
FinalFilteredOrphanSetFileName = "/home/utkinai2/Project1/Orphans_finalfiltered.txt"
FilteredSpacerHitsOrphansFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/orphan_filteredSpacerHits.hits"
FilteredSpacerHitsFullOrphansFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/orphanfull_filteredSpacerHits.hits"
FilteredSpacerHitsFull1OrphansFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/orphanfull_1_filteredSpacerHits.hits"
FilteredSpacerHitsOrphanProteinsFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/orphans_proteins_filteredSpacerHits.hits"
FilteredSpacerHitsBadOrphansFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/bad_orphans_filteredSpacerHits_3-8.hits"


AllOrphansFilteredSpacerHitsFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/Orphans_filtered_SpacerHits.hits"
# SpeciesCountedFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/Species_inOrphanSpacerHits_counted_1.txt"
# FamiliesCountedFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/Families_inOrphanSpacerHits_counted_1.txt"
# ClassesCountedFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/Classes_inOrphanSpacerHits_counted_1.txt"
# OrdersCountedFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/Orders_inOrphanSpacerHits_counted_1.txt"
# PhylaCountedFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/Phyla_inOrphanSpacerHits_counted_1.txt"
# KingdomCountedFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/Kingdoms_inOrphanSpacerHits_counted_1.txt"
SpeciesCountedFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/Species_inBadOrphans_SpacerHits_counted.txt"
FamiliesCountedFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/Families_inBadOrphans_SpacerHits_counted.txt"
ClassesCountedFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/Classes_inBadOrphans_SpacerHits_counted.txt"
OrdersCountedFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/Orders_inBadOrphans_SpacerHits_counted.txt"
PhylaCountedFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/Phyla_inBadOrphans_SpacerHits_counted.txt"
#KingdomCountedFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/Kingdoms_inOrphanSpacerHits_counted_1.txt"
PhylumFileName = "/panfs/pan1/prokdata/db_tmp/all1603.tax.tab"
OrphansWithProteinHitsFileName = "/home/utkinai2/Project1/Orphans_withProtHits.txt"
OrphansWithProteinHitsSpacersFileName = "/home/utkinai2/Project1/Orphans_withProtHits_spacers.txt"
OrphansProteinsSpacerHitsFileName = "/home/utkinai2/Project1/Orphan_ProteinHits_Spacers.hits"
BadOrphansSpacersFileName = "/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/Orphans_with_sameSpacers_spacers.txt"

OrphansWithSameSpacersFileName = "/home/utkinai2/Project1/Orphans_with_sameSpacers.txt"

# Dictionary for ID of filtered Orphan set (final version)
FilteredOrphanID = []
for Line in open(FinalFilteredOrphanSetFileName, "r"):
    LineValues = Line[:-1].split(' ')
    FilteredOrphanID.append(str(LineValues[2].replace('"','') + '_' + LineValues[3] + '_' + LineValues[4]))

print('first block done')

# Dictionaries for  spacers from final filtered sets of partially and fully covered orphans (2 different, because blast of spacers against prokaryotes and viruses has been done separately for them),
# as a key - full spacer ID (f.e. CP000001_292_955_2_spacer_191919_43) as it is shown in blast output, as corresponding value - length of this spacer (in order to
# calculate coverage further)
FullSpacersByID = {}
PartSpacersByID = {}
ProteinSpacersByID = {}
BadOrphanSpacersByID = {}
# with open(GoodFullCoveredSpacersFileName,"r") as f:
#     for line1, line2 in itertools.zip_longest(*[f]* 2):
#         FullSpacersByID[line1[1:-1]] = len(line2)
#
# with open(PartCoveredSpacersFileName, "r") as f:
#     for line1, line2 in itertools.zip_longest(*[f] * 2):
#         PartSpacersByID[line1[1:-1]] = len(line2)

# with open(OrphansWithProteinHitsSpacersFileName, "r") as f:
#     for line1, line2 in itertools.zip_longest(*[f]* 2):
#         ProteinSpacersByID[line1[1:-1]] = len(line2)
#
# with open(BadOrphansSpacersFileName, "r") as f:
#     for line1, line2 in itertools.zip_longest(*[f]* 2):
#         BadOrphanSpacersByID[line1[1:-1]] = len(line2)
print('dictionaries are created')

count = 0
# parsing of blast output and filtering it by the cut-off 95% identity and 95% coverage
# run once, than comment (this part is running during several hours)
#with open(FilteredSpacerHitsBadOrphansFileName, "w") as FilteredSpacerHitsBadOrphansFile:
#with open(FilteredSpacerHitsFull1OrphansFileName, "w") as FilteredSpacerHitsFull1OrphansFile:
with open(FilteredSpacerHitsOrphanProteinsFileName, "w") as FilteredSpacerHitsOrphanProteinsFile:
    for line in open("/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/Orphans_withProtHits_spacers.txt_vs_all1603.hits" ,"r"):
        if line[0] == "#":
            continue
        LineValues = line[:-1].split('\t')
        ID = re.split(r'_[0-9]{1,3}_spacer', LineValues[0])
        matches = sum(nt1 == nt2 for nt1, nt2 in zip(LineValues[6], LineValues[7]))
        identity = 100.0*matches/len(LineValues[6])
        coverage = 100.0*len(LineValues[6])/ProteinSpacersByID[LineValues[0]]
        if identity >= 95.0 and coverage >= 95.0:
            print('YES' , ID[0], ' ', identity, '% ', coverage, '%')
            FilteredSpacerHitsOrphanProteinsFile.write(line)
    # for number in range(1,11):
    #     for line in io.open("/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/orphan" + str(number) + "_vs_all1603.hits" , "r"):
    #         print('FILE', number)
    #         if line[0] == "#":
    #             continue
    #         LineValues = line[:-1].split('\t')
    #         ID = re.split(r'_[0-9]{1,3}_spacer',LineValues[0])
    #         if ID[0] in FilteredOrphanID:
    #             matches = sum(nt1 == nt2 for nt1, nt2 in zip(LineValues[6], LineValues[7]))
    #             identity = 100.0*matches/len(LineValues[6])
    #             coverage = 100.0*len(LineValues[6])/PartSpacersByID[LineValues[0]]
    #             if identity >= 95.0 and coverage >= 95.0:
    #                 print('YES' , ID[0], ' ', identity, '% ', coverage, '%')
    #                 FilteredSpacerHitsOrphansFile.write(line)
    #for number in range(1, 9):
     # for number in range(5,7):
     #    for line in io.open("/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/orphanfull" + str(number) + "_vs_all1603.hits", "r"):
     #        print('FILE', number)
     #        if line[0] == "#":
     #            continue
     #        LineValues = line[:-1].split('\t')
     #        ID = re.split(r'_[0-9]{1,3}_spacer', LineValues[0])
     #        if ID[0] in FilteredOrphanID:
     #            matches = sum(nt1 == nt2 for nt1, nt2 in zip(LineValues[6], LineValues[7]))
     #            identity = 100.0 * matches / len(LineValues[6])
     #            coverage = 100.0 * len(LineValues[6]) / FullSpacersByID[LineValues[0]]
     #            if identity >= 95.0 and coverage >= 95.0:
     #                print('YES' , ID[0], ' ', identity, '% ', coverage, '%')
     #                FilteredSpacerHitsFull1OrphansFile.write(line)
    # for number in range(7,9):
    #    for line in io.open("/panfs/pan1/orphancrispr/SpacerHits_against_prokandvira/bad_orphans" + str(number) + "_vs_all1603.hits", "r"):
    #        print('FILE', number)
    #        if line[0] == "#":
    #            continue
    #        LineValues = line[:-1].split('\t')
    #        ID = re.split(r'_[0-9]{1,3}_spacer', LineValues[0])
    #        matches = sum(nt1 == nt2 for nt1, nt2 in zip(LineValues[6], LineValues[7]))
    #        identity = 100.0 * matches / len(LineValues[6])
    #        coverage = 100.0 * len(LineValues[6]) / BadOrphanSpacersByID[LineValues[0]]
    #        if identity >= 95.0 and coverage >= 95.0:
    #            print('YES' , ID[0], ' ', identity, '% ', coverage, '%')
    #            FilteredSpacerHitsBadOrphansFile.write(line)

# Creating a dictionary, all species available from 1603 database as keys and corresponding list [family, class, order, phylum, kingdom] as value
SpeciesGeneralDict = {}
for line in open(PhylumFileName, "r"):
    LineValues = line[:-1].split('\t')
    if LineValues[-1] not in SpeciesGeneralDict:
        SpeciesGeneralDict[re.sub(r'\\[|\\]', "",LineValues[-1])] = [LineValues[-3], LineValues[-4], LineValues[-5], LineValues[-6], LineValues[-7]]


# After blastn agains entire prokaryotic db, spacer hits from orphan arrays filtered by cut-off 95% identity and 95% coverage are counted by species,
# and creating a dictionary, keys - organism ID (LK1111_11_111), values - origins of spacers (species) within this array
SpeciesSpacerHitsInArray = {}
LastSpacer = 12
SpeciesOneSpacer = []


#for line in open(AllOrphansFilteredSpacerHitsFileName, "r"):
for line in open(FilteredSpacerHitsBadOrphansFileName, "r"):
    LineValues = line[:-1].split('\t')
    SpacerID = LineValues[0].split('_')
    Array = str(SpacerID[0]+'_'+SpacerID[1]+'_'+SpacerID[2])
    if LineValues[14].split(' ')[1] == "sp.":
        Specie = re.sub(r'\\[|\\]', "", str(LineValues[14].split(' ')[0] + ' ' + LineValues[14].split(' ')[1] + ' ' + LineValues[14].split(' ')[2]).strip(','))
    else:
        Specie = re.sub(r'\\[|\\]', "", str(LineValues[14].split(' ')[0] + ' ' + LineValues[14].split(' ')[1]).strip(','))
    if not re.match(r'Length', Specie):
        if Array in SpeciesSpacerHitsInArray:
            if Specie not in SpeciesSpacerHitsInArray[Array]:
                SpeciesSpacerHitsInArray[Array] += [Specie]
        else:
            SpeciesSpacerHitsInArray[Array] = [Specie]
print(len(SpeciesSpacerHitsInArray.keys()))

FamiliesSpacerHitsInArray = {}
ClassesSpacerHitsInArray = {}
OrdersSpacerHitsInArray = {}
PhylaSpacerHitsInArray = {}

count1 = 0
for Array in SpeciesSpacerHitsInArray.keys():
    count += 1
    try:
        for Specie in SpeciesSpacerHitsInArray[Array]:
            for key in SpeciesGeneralDict.keys():
                # if count >294:
                #     print(key)
                if key != '': #and (Specie!= '' or Specie != ']'):
                    if re.match(Specie,key):
                        print(count)
                        print(Specie, ', ', key)

        # if (len(key.split(' ')) == 2 and Specie == str(key.split(' ')[0] + ' ' + key.split(' ')[1])) or \
        #         (len(key.split(' ')) > 2 and Specie == str(key.split(' ')[0] + ' ' + key.split(' ')[1] + ' ' + key.split(' ')[2])) :
                        if Array not in FamiliesSpacerHitsInArray:
                            FamiliesSpacerHitsInArray[Array] = [SpeciesGeneralDict[key][0]]
                        else:
                            if SpeciesGeneralDict[key][0] not in FamiliesSpacerHitsInArray[Array]:
                                FamiliesSpacerHitsInArray[Array] += [SpeciesGeneralDict[key][0]]
                        if Array not in ClassesSpacerHitsInArray:
                            ClassesSpacerHitsInArray[Array] = [SpeciesGeneralDict[key][1]]
                        else:
                            if SpeciesGeneralDict[key][1] not in ClassesSpacerHitsInArray[Array]:
                                ClassesSpacerHitsInArray[Array] += [SpeciesGeneralDict[key][1]]
                        if Array not in OrdersSpacerHitsInArray:
                            OrdersSpacerHitsInArray[Array] = [SpeciesGeneralDict[key][2]]
                        else:
                            if SpeciesGeneralDict[key][2] not in OrdersSpacerHitsInArray[Array]:
                                OrdersSpacerHitsInArray[Array] += [SpeciesGeneralDict[key][2]]
                        if Array not in PhylaSpacerHitsInArray:
                            PhylaSpacerHitsInArray[Array] = [SpeciesGeneralDict[key][3]]
                        else:
                            if SpeciesGeneralDict[key][3] not in PhylaSpacerHitsInArray[Array]:
                                PhylaSpacerHitsInArray[Array] += [SpeciesGeneralDict[key][3]]
    except:
        print(key)
#print(PhylaSpacerHitsInArray)

with open(PhylaCountedFileName, "w") as PhylaCountedFile:
    with open(OrdersCountedFileName, "w") as OrdersCountedFile:
        with open(ClassesCountedFileName, "w") as ClassesCountedFile:
            with open(FamiliesCountedFileName, "w") as FamiliesCountedFile:
                with open(SpeciesCountedFileName, "w") as SpeciesCountedFile:
                    for Array in SpeciesSpacerHitsInArray:
                        SpeciesCountedFile.write(Array.split('_')[0] + '\t' + Array.split('_')[1] + '\t' + Array.split('_')[2] + '\t' + str(SpeciesSpacerHitsInArray[Array]) + '\t' + str(len(SpeciesSpacerHitsInArray[Array])) + '\n')
                    for Array in FamiliesSpacerHitsInArray:
                        FamiliesCountedFile.write(Array.split('_')[0] + '\t' + Array.split('_')[1] + '\t' + Array.split('_')[2] + '\t' + str(FamiliesSpacerHitsInArray[Array]) + '\t' + str(len(FamiliesSpacerHitsInArray[Array])) + '\n')
                    for Array in ClassesSpacerHitsInArray:
                        ClassesCountedFile.write(Array.split('_')[0] + '\t' + Array.split('_')[1] + '\t' + Array.split('_')[2] + '\t' + str(ClassesSpacerHitsInArray[Array]) + '\t' + str(len(ClassesSpacerHitsInArray[Array])) + '\n')
                    for Array in OrdersSpacerHitsInArray:
                        OrdersCountedFile.write(Array.split('_')[0] + '\t' + Array.split('_')[1] + '\t' + Array.split('_')[2] + '\t' + str(OrdersSpacerHitsInArray[Array]) + '\t' + str(len(OrdersSpacerHitsInArray[Array])) + '\n')
                    for Array in PhylaSpacerHitsInArray:
                        PhylaCountedFile.write(Array.split('_')[0] + '\t' + Array.split('_')[1] + '\t' + Array.split('_')[2] + '\t' + str(PhylaSpacerHitsInArray[Array]) + '\t' + str(len(PhylaSpacerHitsInArray[Array])) + '\n')


# print(len(FamilyDict.keys()))
# print(len(ClassDict.keys()))
# print(len(OrderDict.keys()))
# print(len(PhylumDict.keys()))
