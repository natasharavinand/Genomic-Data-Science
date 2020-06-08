#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 10:10:44 2020

@author: natasharavinand
"""

''' This file contains 3 naive matching algorithms for fittings reads onto a DNA sequence. It includes:
    
    - a basic naive matching algorithm (naive())
    - a naive algorithm that allows up to 2 mismatches (naive_2mm())
    - a strand-aware naive algorithm that additionally matches reads according to their reverse complement 
    (strand_aware_naive())

Additionally, this file includes programs to:
    
    - Read and parse a fastq file (read_fastq())
    - Find the reverse complement of a sequence (reverse_complement())
    - Find where N occurs most in a file of reads to determine an inaccurate sequencing cycle (find_n_by_pos())
    - Find the longest common prefix of two sequences (longest_common_prefix())
    - Turn Q into Phred++33 ASCII-encoded quality (QtoPhred33())
    - Turn Phred+33 ASCII-encoded quality into Q (phred33ToQ())
    - Generate a list of Q frequencies from Phred++33 to plot as a histogram (create_list_for_hist())
    - Find GC content by position of reads (find_gc_by_pos())
    
'''
    

def naive(p, r):
    ''' This is a basic naive matching algorithm that takes in a pattern p and compares it to reference r '''
    occurrences = []
    for i in range(len(r) - len(p) + 1): # Loop over alignments
        match = True
        for j in range(len(p)): # Look over characters
            if r[i+j] != p[j]: # Compare characters
                match = False # Record mismatch
                break
        if match:
            occurrences.append(i) # Record match
    return occurrences



def naive_2mm(p, r):
    ''' This is a basic naive matching algorithm that allows 2 mismatches and
    takes in a pattern p and compares it to reference r '''
    occurrences = []
    for i in range(len(r) - len(p) + 1): 
        mismatches = 0
        match = True
        for j in range(len(p)): 
            if r[i+j] != p[j]: 
                if mismatches > 1: #if more than two mismatches
                    match = False 
                    break
                else:
                    mismatches += 1
            else:
                continue
        if match:
            occurrences.append(i) # Record match
    return occurrences



def strand_aware_naive(p, r):
    ''' This is a strand-aware naive algorithm that additionally matches reads according to their reverse complement using
    pattern p and reference r '''
    rc = reverse_complement(p) #reverse complement of p
    occurrences = []
    
    for i in range(len(r) - len(p) + 1):
        match = True
        for j in range(len(p)):
            if r[i+j] != p[j]:
                match = False
                break
        if match:
            occurrences.append(i)
            
    for i in range(len(r) - len(rc) + 1):
        r_match = True
        for j in range(len(rc)):
            if r[i+j] != rc[j]:
                r_match = False
                break
        if r_match and p != rc:
            occurrences.append(i)
            
    return occurrences



def read_fastq(filename):
    ''' This function reads and parses a fastq file and returns a list of sequences and base qualities '''
    sequences = []
    qualities = []
    
    with open(filename) as f:
        while True:
            f.readline()
            seq = f.readline().rstrip()
            f.readline()
            qual = f.readline().rstrip()
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
            
    return sequences, qualities



def reverse_complement(s):
    ''' This function returns the reverse complement of sequence s '''
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    t = ''
    for base in s:
        t = complement[base] + t #this causes us to go from the back to the front
    return t



def find_n_by_pos(reads):
    ''' This function finds where N occurs most in a file of reads to determine an inaccurate 
    sequencing cycle '''
    n_pos = []
    found = 0
    
    for read in reads:
        for i in range(len(read)):
            if read[i] == "N":
                found += 1
                n_index = i
                n_pos.append(n_index)
                
    return max(set(n_pos), key=n_pos.count)
                


def longest_common_prefix(s1, s2):
    ''' This function returns the longest common prefix of two sequences s1 and s2 '''
    i = 0
    while i < len(s1) and i < len(s2) and s1[i] == s2[i]:
        i += 1
    return s1[:i]



def QtoPhred33(Q):
    """ Turn Q into Phred++33 ASCII-encoded quality """
    return chr(Q + 33)



def phred33ToQ(qual):
    """ Turn Phred+33 ASCII-encoded quality into Q """
    return ord(qual) - 33



def create_list_for_hist(qualities):
    ''' This function generates a list of Q frequencies from Phred++33 to plot as a histogram '''
    hist = [0] * len(qualities)
    for qual in qualities:
        for phred in qual:
            q = phred33ToQ(phred)
            hist[q] += 1
    return hist



def find_gc_by_pos(reads):
    ''' This function finds GC content by position of reads given '''
    gc = [0] * len(reads)
    totals = [0] * len(reads)
    
    for read in reads:
        for i in range(len(read)):
            if read[i] == 'C' or read[i] == 'G':
                gc[i] += 1
            totals[i] += 1
            
    for i in range(len(gc)):
        if totals[i] > 0:
            gc[i] /= float(totals[i])
            
    return gc
