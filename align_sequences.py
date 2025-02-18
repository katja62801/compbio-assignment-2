#!/usr/bin/env python

"""
    usage:
        align_sequences [options] seq1.fa seq2.fa
    where the options are:
        -h,--help : print usage and quit
        -m,--match: score of a match in the alignment [2]
        -x,--mismatch: penalty for a mismatch in the alignment [1]
        -g,--gapopen: penalty for opening a new gap [4]
        -e,--gapextend: penalty for extending a gap [1]
"""

from sys import argv, stderr
from getopt import getopt, GetoptError

# a simple function to read the name and sequence from a file
# The file is expected to have just one contig/sequence. This function
# checks the assumption and complains if it is not the case.
def read_single_contig_fasta(filename):
    names = []
    sequences = []
    with open(filename, 'r') as f:
        line = f.readline()
        assert line.startswith(">")
        names.append(line.strip().split("\t"))
        sequence = ""
        for line in f:
            if line.startswith(">"):
                sequences.append(sequence)
                names.append(line.strip().split("\t"))
                sequence = ""
            else:
                for x in line.strip():
                    if x not in ["A", "C", "G", "T"]:
                        print("Unknown nucleotide {}".format(x), file=stderr)
                        exit(3)
                sequence += line.strip()

    sequences.append(sequence)
    assert len(names) == 1
    assert len(sequences) == 1
    return names[0], sequences[0]

def smith_waterman(seq1, seq2, match, mismatch, gapopen, gapextend):

    import numpy as np

    #create matrix full of zeros to start
    m = len(seq1)
    n = len(seq2)
    matrix = numpy.zeros(m+1,n+1)
    positionmatrix = numpy.zeros(m+1,n+1)
    max_score = 0
    maxposition = (i,j)

    for i in range (2, m+1):
        for j in range (2, n+1):
        
            if seq1[i-1] == seq2[j-1]:
                diagonalscore = matrix[i-1,j-1] + match
            else:
                diagonalscore = matrix[i-1, j-1] - mismatch
        
            if seq1opengap = true:
                seq1gapscore = matrix[i-1, j] - gapextend
            else: 
                seq1gapscore = matrix[i-1, j] - gapopen
        
            if seq2opengap = true: 
                seq2gapscore = matrix[i, j-1] - gapextend
            else: 
                seq2gapscore = matrix[i, j-1] - gapopen
        
            matrix[i,j] = max(diagonalscore, seq1gapscore, seq2gapscore)
            if matrix[i,j] > maxscore:
                max_score = matrix[i,j]
                maxposition = (i,j)

            if matrix[i,j] = seq1gapscore:
                seq1opengap = true
                positionmatrix[i,j] = 1 #1 means there is a gap in seq1
            elif matrix[i,j] = seq2gapscore:
                seq2opengap = true
                positionmatrix[i,j] = 2 #2 means there is a gap in seq2
            else:
                positionmatrix[i,j] = 4 #4 means there is no gap
    
    alnseq1 = ""
    alnseq2 = ""

    i, j = maxposition
    while i > 0 and j > 0 and matrix[i, j] > 0:
        if positionmatrix[i,j] = 4:
            alnseq1 = seq1[i-1] + alnseq1
            alnseq2 = seq2[j-1] + alnseq2
            i = i-1
            j = j-1
        elif positionmatrix[i,j] = 2:
            alnseq1 = seq1[i-1] + alnseq1
            alnseq2 = "-"
            i = i-1
        elif positionmatrix[i,j] = 1:
            alnseq1 = "-" + alnseq1
            alnseq2 = seq2[j-1] + alnseq2
            j = j-1
            
    


    return max_score, alnseq1, alnseq2
    

def main(filename1, filename2, match, mismatch, gapopen, gapextend):
    # read the name and sequence from the file
    name1, seq1 = read_single_contig_fasta(filename1)
    name2, seq2 = read_single_contig_fasta(filename2)

    # this function takes as input two nucleotide sequences along with
    # scores for an alignment match, mismatch, opening a new gap, and 
    # extending an existing gap. This should return the maximum alignment
    # score as well as the alignment. For examples see the testdriver script
    max_score, alnseq1, alnseq2 = smith_waterman(seq1, seq2, 
                                  match, mismatch, gapopen, gapextend)
    
    print("Maximum alignment score: {}".format(max_score))
    print("Sequence1 : {}".format(alnseq1))
    print("Sequence2 : {}".format(alnseq2))

if __name__ == "__main__":
    try:
        opts, args = getopt(argv[1:],
                     "hm:x:g:e:",
                     ["help", "match=", "mismatch=", "gapopen=", "gapextend="])
    except GetoptError as err:
        print(err)
        print(__doc__, file=stderr)
        exit(1) 

    match = 2
    mismatch = 1
    gapopen = 4
    gapextend = 1

    for o, a in opts:
        if o in ("-h", "--help"):
            print(__doc__, file=stderr)
            exit()
        elif o in ("-m", "--match"):
            match = float(a)
        elif o in ("-x", "--mismatch"):
            mismatch = float(a)
        elif o in ("-g", "--gapopen"):
            gapopen = float(a)
        elif o in ("-e", "--gapextend"):
            gapextend = float(a)
        else:
            assert False, "unhandled option"

    if len(args) != 2:
        print(__doc__, file=stderr)
        exit(2)

    main(args[0], args[1], match, mismatch, gapopen, gapextend)
