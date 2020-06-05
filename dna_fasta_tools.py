#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 13:41:41 2020

@author: natasharavinand
"""
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq

class dna_fasta_tools():
    
    '''This class is created to serve as a set of tools to utilize in simple analysis of DNA. It works with
    inputs of DNA sequences in multi-FASTA format. 
    
    This class can perform analysis on the following operations:
        
        - Identifying the number of records in a file (num_records())
        -Aligning desired sequence from sequence identifier to NCBI BLAST and returning the results under a provided e-value threshold
        -Identifying the lengths of sequences from a sequence dictionary (length_sequences())
        -Identifying open reading frames (ORFs) from a given reading frame (identify_orfs())
        -Finding repeats composed of all possible samples of n-length and their 
        respective counts (sequence_repeats())
        -Finding sequences of repeats from an entire file and their n-length counts (repeat_finder())
        -Generate a corresponding protein sequence from a sequence provided by a sequence identifier
    '''
    
    
    def __init__(self, file):
        ''' Initializes an object of the dna_fasta_tools class, such that each object has its file attached to it
        as well as a dictionary that is critical in usage for functions below'''
        self.dict = {}
        self.file = file
        try:
            reader = open(self.file)
        except IOError:
            raise IOError("File was not found.\n")
        for line in reader:
            #discard newline at the end
            line = line.rstrip()
            #distinguish header from sequence
            if line.startswith('>'): #if header
                words = line.split()
                name = words[0][1:]
                self.dict[name] = ""
            else: #if sequence
               self.dict[name] += line
        reader.close()
        
        
        
    def num_records(self):
        '''Identify the number of records in this file from a given sequence dictionary. Prints out
        the number'''
        num_records = len(self.dict)
        print("The number of records in this file is:", num_records, "\n")
        
        
        
    def seq_id_blast_alignment(self, seq_id, e_val):
        '''Takes in a sequence identifier (name of sequence) and e-value threshold. Compares it to the NCBI BLAST alignment tool.
        Returns values found under given e_value threshold. Request goes through NCBIWWW server so may take
        a few minutes.'''
        if type(seq_id) != str:
            raise ValueError("Please input the sequence identifier as a string.")
        seq = self.dict[seq_id]
        sequence = Seq(seq)
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
        blast_record = NCBIXML.read(result_handle)
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < e_val:
                    print("***ALIGNMENT***")
                    print("sequence: ", alignment.title)
                    print("length: ", alignment.length)
                    print("e value: ", hsp.expect)
                    print(hsp.query)
                    print(hsp.match)
                    print(hsp.sbjct)

        
        
    def length_sequences(self):
        '''Identify the lengths of sequences from a given sequence dictionary. Prints out the
        longest sequence with x number of records as well as the shortest sequence with x number
        of records '''
        lengths = {}
        for seq in self.dict:
            lengths[seq] = len(self.dict[seq])
        #Sorting dictionary from greatest to least in terms of lengths
        ordered_lengths = []
        for k in sorted(lengths, key=lengths.get, reverse=True):
            ordered_lengths.append([k, lengths[k]])
        low_counter = 1
        high_counter = 1
        index = 0
        #checking to see if there are multiple longest lengths
        while True:
            if ordered_lengths[index][1] == ordered_lengths[index + 1][1] and index < \
                len(ordered_lengths):
                    high_counter += 1
            else:
                break
        #checking to see if there are multiple shortest lengths
        index = len(ordered_lengths) - 1
        while True:
            if ordered_lengths[index][1] == ordered_lengths[index-1][1] and index > 0:
                low_counter += 1
            else:
                break
        print("The longest sequence is ", ordered_lengths[0][0], "with", high_counter, "record(s) of length",\
              ordered_lengths[0][1], "and the shortest sequence is ", ordered_lengths[len(ordered_lengths)-1][0], "with", \
                  low_counter, "record(s) of length", ordered_lengths[len(ordered_lengths) - 1][1], ".", "\n")
            
            
            
    def identify_orfs(self, frame, sort=False, seq_id=""):
        '''Identify ORFs (open reading frames) in this file given a reading frame (1, 2, or 3)
        and a sequence dictionary. Returns a dictionary of ORFs with a key of the ORF string
        and the sequence identifier, starting position, length of the ORF, and frame used as a 
        list value. Works only in 5' to 3' direction
        
        "frame" –> reading frame (either 1,2,3) to read the frame in 
        "sort" –> if true, will print the greatest to least length in ORF, its sequence identifier, and its starting position
        "seq_id" –> If sequence id given, will print the longest ORF with a given sequence identifier, its length, and position
        '''
        stop_codons = ['TAA', 'TAG', 'TGA']
        orfs = {}
        if frame != 1 and frame != 2 and frame != 3:
            raise ValueError("Please input either 1, 2, or 3 as a reading frame.")
        for seq in self.dict:
            sequence = self.dict[seq]
            for i in range(frame-1, len(sequence), 3):
                codon = sequence[i:i+3]
                if codon == 'ATG': #if start codon
                    starting_pos = i + 1
                    is_orf = False
                    potential_orf = "ATG"
                    for j in range(i+3, len(sequence), 3):
                        next_codon = sequence[j:j+3]
                        if next_codon in stop_codons: #a stop codon
                            potential_orf += next_codon
                            is_orf = True
                            break
                        else:
                            potential_orf += next_codon
                    if is_orf:
                        length = len(potential_orf)
                        orfs[potential_orf] = [seq, starting_pos, length, frame]
                else:
                    continue
        if sort == True:
            sorted_list = sorted(orfs.items(), key=lambda x:x[1][2], reverse=True)
            longest_orf = next(iter(sorted_list))
            print("The longest ORF in reading frame,", longest_orf[1][3], "is in sequence identifier ", \
                  longest_orf[1][0], "with length", longest_orf[1][2], "in position", longest_orf[1][1], ".\n")
        if seq_id != "":
            sorted_list = sorted(orfs.items(), key=lambda x:x[1][2], reverse=True)
            for val in sorted_list:
                if val[1][0] == seq_id:
                    print("The longest ORF for sequence identifier", val[1][0], "has length", \
                          val[1][2], "starting at position", val[1][1],".\n") 
                    break
        return orfs
        
    
    
    def sequence_repeats(n, sequence):
        '''Takes in a sequences and creates a repeat dictionary composed of all possible samples of n-length
        and their respective counts'''
        repeats_dict = {}
        cap = n
        for i in range(0, len(sequence)):
            repeat = sequence[i:i+n]
            if cap == len(repeat):
                if repeat in repeats_dict:
                    repeats_dict[repeat] += 1
                else:
                    repeats_dict[repeat] = 1
        return repeats_dict
    
    
    
    def repeat_finder(self, n):
        '''Takes in a n-length and uses helper function sequence_repeats() to return a dictionary of total
        repeats of n-length across the file'''
        seq_repeats = {}
        for name, sequence in self.dict.items():
            repeats_dict = dna_fasta_tools.sequence_repeats(n, sequence)
            seq_repeats[name] = repeats_dict
        total_repeats = {}
        for dictionary in seq_repeats.values():
            for name in dictionary:
                if name in total_repeats:
                    total_repeats[name] += dictionary[name]
                else:
                    total_repeats[name] = dictionary[name]
        return total_repeats
    
    
    
    def to_protein(self, seq_id):
        '''Translates a given DNA sequence with a sequence identifier into a protein'''
        seq = self.dict[seq_id]
        codon_table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
        }
        protein_sequence = ""
        if len(seq) % 3 != 0:
            raise ValueError("Please enter a DNA sequence with complete codons of length 3.")
        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            protein_element = codon_table[codon]
            protein_sequence += protein_element
        return protein_sequence
        

        