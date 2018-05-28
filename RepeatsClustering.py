# import numpy
# import h5py
# import psutil
# import natsort
# import prody
# import pyRMSD
# import bottleneck

# from cogent import LoadSeqs
# from cogent.phylo import distance
# from cogent.evolve.models import HKY85

import itertools
import sklearn
import re
import csv
import pandas as pd
import numpy
from sklearn.cluster import AffinityPropagation
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from pandas import DataFrame
from Bio.SubsMat import MatrixInfo as matlist
from sklearn.metrics.pairwise import pairwise_distances

def zeros(shape):
    retval = []
    for x in range(shape[0]):
        retval.append([])
        for y in range(shape[1]):
            retval[-1].append(0)
    return retval


match_award = 10
mismatch_penalty = -3
gap_penalty = -5  # both for opening and extanding


def match_score(alpha, beta):
    if alpha == beta:
        return match_award
    elif alpha == '-' or beta == '-':
        return gap_penalty
    else:
        return mismatch_penalty


def finalize(align1, align2):
    align1 = align1[::-1]  # reverse sequence 1
    align2 = align2[::-1]  # reverse sequence 2

    i, j = 0, 0

    # calcuate identity, score and aligned sequeces
    symbol = ''
    found = 0
    score = 0
    identity = 0
    for i in range(0, len(align1)):
        # if two AAs are the same, then output the letter
        if align1[i] == align2[i]:
            symbol = symbol + align1[i]
            identity = identity + 1
            score += match_score(align1[i], align2[i])

        # if they are not identical and none of them is gap
        elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-':
            score += match_score(align1[i], align2[i])
            symbol += ' '
            found = 0

        # if one of them is a gap, output a space
        elif align1[i] == '-' or align2[i] == '-':
            symbol += ' '
            score += gap_penalty

    identity = float(identity) / len(align1) * 100
    return(score)
    # print('Identity =', "%3.3f" % identity, 'percent')
    print ('Score =', score)
    # print(align1)
    # print (symbol)
    # print(align2)


def needle(seq1, seq2):
    m, n = len(seq1), len(seq2)  # length of two sequences

    # Generate DP table and traceback path pointer matrix
    score = zeros((m + 1, n + 1))  # the DP table

    # Calculate DP table
    for i in range(0, m + 1):
        score[i][0] = gap_penalty * i
    for j in range(0, n + 1):
        score[0][j] = gap_penalty * j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + match_score(seq1[i - 1], seq2[j - 1])
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            score[i][j] = max(match, delete, insert)

    # Traceback and compute the alignment
    align1, align2 = '', ''
    i, j = m, n  # start from the bottom right cell
    while i > 0 and j > 0:  # end toching the top or the left edge
        score_current = score[i][j]
        score_diagonal = score[i - 1][j - 1]
        score_up = score[i][j - 1]
        score_left = score[i - 1][j]

        if score_current == score_diagonal + match_score(seq1[i - 1], seq2[j - 1]):
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif score_current == score_left + gap_penalty:
            align1 += seq1[i - 1]
            align2 += '-'
            i -= 1
        elif score_current == score_up + gap_penalty:
            align1 += '-'
            align2 += seq2[j - 1]
            j -= 1

    # Finish tracing up to the top left cell
    while i > 0:
        align1 += seq1[i - 1]
        align2 += '-'
        i -= 1
    while j > 0:
        align1 += '-'
        align2 += seq2[j - 1]
        j -= 1

    return finalize(align1, align2)


def water(seq1, seq2):
    m, n = len(seq1), len(seq2)  # length of two sequences

    # Generate DP table and traceback path pointer matrix
    score = zeros((m + 1, n + 1))  # the DP table
    pointer = zeros((m + 1, n + 1))  # to store the traceback path

    max_score = 0  # initial maximum score in DP table
    # Calculate DP table and mark pointers
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score_diagonal = score[i - 1][j - 1] + match_score(seq1[i - 1], seq2[j - 1])
            score_up = score[i][j - 1] + gap_penalty
            score_left = score[i - 1][j] + gap_penalty
            score[i][j] = max(0, score_left, score_up, score_diagonal)
            if score[i][j] == 0:
                pointer[i][j] = 0  # 0 means end of the path
            if score[i][j] == score_left:
                pointer[i][j] = 1  # 1 means trace up
            if score[i][j] == score_up:
                pointer[i][j] = 2  # 2 means trace left
            if score[i][j] == score_diagonal:
                pointer[i][j] = 3  # 3 means trace diagonal
            if score[i][j] >= max_score:
                max_i = i
                max_j = j
                max_score = score[i][j]

    align1, align2 = '', ''  # initial sequences

    i, j = max_i, max_j  # indices of path starting point

    # traceback, follow pointers
    while pointer[i][j] != 0:
        if pointer[i][j] == 3:
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif pointer[i][j] == 2:
            align1 += '-'
            align2 += seq2[j - 1]
            j -= 1
        elif pointer[i][j] == 1:
            align1 += seq1[i - 1]
            align2 += '-'
            i -= 1
    finalize(align1, align2)


ScoreMatrixFileName = "/home/utkinai2/Project1/repeats_scorematrix.csv"
RepeatsFastaFileName = "/home/utkinai2/Project1/ClusteringRepeats/identified_repeats.txt"
count = 0

# data = pd.DataFrame({'repeat_seq' : [], 'ID': []})
RepeatsDict = {'ID':[], 'repeat_seq' : []}
with open(RepeatsFastaFileName, "r") as f:
    # for line1, line2 in itertools.zip_longest(*[f]* 2):
    for line in f:
        LineValues = line[:-1].split(' ')
        Match = re.search(r'[0-9]{1,}', LineValues[0])
        if Match:
            RepeatsDict['ID'].append(Match.group(0))
        else:
            RepeatsDict['repeat_seq'].append(line)
RepeatsTable = DataFrame.from_dict(RepeatsDict, orient='columns', dtype=None)
RepeatsTable['repeat_seq'].replace(regex=True,inplace=True,to_replace=r'\n',value=r'')
# RepeatsTable.index.name = 'ID'
# RepeatsTable.reset_index()
# print(RepeatsTable)

ScoreMatrix = pd.DataFrame(index = RepeatsDict['ID'] , columns = RepeatsDict['ID'])
# print(ScoreMatrix.shape)
SimilarityMatrix = numpy.zeros(shape=(len(RepeatsDict['ID']),len(RepeatsDict['ID'])))


matrix = matlist.blosum62
Results = []
#with open(ScoreMatrixFileName, "w") as ScoreMatrixFile:
for seq0, seq1 in itertools.combinations(RepeatsTable['repeat_seq'], 2):
    AlignmentScore = pairwise2.align.globalms(seq0,seq1,1,-0.1,-1,-0.1) #, score_only=True)
    Results.append(needle(seq0,seq1))

    # if str(AlignmentScore).isdigit() == False:
    #     print(RepeatsTable.loc[RepeatsTable['repeat_seq'] == seq0, 'ID'].iloc[0], ' ',
    #           RepeatsTable.loc[RepeatsTable['repeat_seq'] == seq1, 'ID'].iloc[0])
    BestHit = AlignmentScore[0]
    Score = BestHit[2]
    # SimilarityMatrix[RepeatsTable.loc[RepeatsTable['repeat_seq'] == seq0, 'ID'].iloc[0], RepeatsTable.loc[RepeatsTable['repeat_seq'] == seq1, 'ID'].iloc[0]] = float(Score)
    # SimilarityMatrix[RepeatsTable.loc[RepeatsTable['repeat_seq'] == seq1, 'ID'].iloc[0], RepeatsTable.loc[RepeatsTable['repeat_seq'] == seq0, 'ID'].iloc[0]] = float(Score)
    ScoreMatrix.at[RepeatsTable.loc[RepeatsTable['repeat_seq'] == seq0, 'ID'].iloc[0],
                   RepeatsTable.loc[RepeatsTable['repeat_seq'] == seq1, 'ID'].iloc[0]] = Score #needle(seq0,seq1)
    ScoreMatrix.at[RepeatsTable.loc[RepeatsTable['repeat_seq'] == seq1, 'ID'].iloc[0],
                   RepeatsTable.loc[RepeatsTable['repeat_seq'] == seq0, 'ID'].iloc[0]] = Score  #needle(seq0,seq1)
    # ScoreMatrixFile.write(ScoreMatrix)
# print(SimilarityMatrix)
#
# # pairwise_distances(RepeatsTable['repeat_seq'])
print(ScoreMatrix)
#print(Results)



