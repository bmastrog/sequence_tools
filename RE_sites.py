#!/usr/bin/python3

"""
This script provides functions that allow for analysis of restriction enzyme sites.  

One function determines nucleotide compositions of DraII sites (RGGNCCY) in a given 
set of DNA sequences, where R is a purine, N is any nucleotide, and Y is a pyrimadine.

The second function finds start and end positions of styI enzyme sites in a sequence 
(CC[AT][AT]GG).

"""

import argparse
import re

if  __name__== '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="input fasta file w. nt sequences")
    parser.add_argument("-e", help="RE of interest")
    args = parser.parse_args() 
    
def load_FASTA(filename):
    '''INPUT: a string representing a file name.  The file can contain one or more fasta files.
       OUTPUT: the sequence or sequences as a string as a list'''

    with open (filename, 'r') as f:
        contents = f.read() #reads the file as a string
        
        allsequences = contents.split('>') 
        sequences = []  #initialize storage variable
        sequences = [seq.partition('\n')[2] for seq in allsequences] 
        #creates strings before \n, \n, and after \n for break/join point
        sequences = [''.join(seq.split('\n')) for seq in sequences]
    
    return sequences

def seq_statistic(sequence):
   '''Takes input sequences and returns sequence composition of DraII sites'''
   DraII = "[AG]GG[ACTG]CC[CT]" 
   sequence = ''.join(sequence) 
   found = re.findall(DraII,sequence) 
   allsites = []   
   for i in found: 
      allsites.append(i) 
      n = len(allsites[0]) 
      A = [0]*n 
      T = [0]*n
      G = [0]*n
      C = [0]*n
   for dna in allsites: 
      for index, base in enumerate(dna):  
         if base == 'A': 
            A[index] += 1 
         elif base == 'C': 
            C[index] += 1 
         elif base == 'G': 
            G[index] += 1 
         elif base == 'T': 
            T[index] += 1       
   Aa = [round(pos/float(len(allsites)),2) for pos in A] 
   Cc = [round(pos/float(len(allsites)),2) for pos in C]
   Gg = [round(pos/float(len(allsites)),2) for pos in G]
   Tt = [round(pos/float(len(allsites)),2) for pos in T]
    
   print('\t1\t2\t3\t4\t5\t6\t7')
   print('A\t', str(Aa).replace(',','\t') )
   print('C\t', str(Cc).replace(',','\t'))
   print('G\t', str(Gg).replace(',','\t'))
   print('T\t', str(Tt).replace(',','\t'))

sequence = load_FASTA(args.f)

if args.e == 'draII':
   seq_statistic(sequence)


def ReverseComplement1(seq):
   "Input sequences and it outputs the reverse complement strand"
   seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':"N"}
   return "".join([seq_dict[base] for base in reversed(seq)])


def restriction_site_scan(sequence):     
    sequence = ''.join(sequence)
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','N':'N'} 
    comp = [complement[i] for i in sequence] 
    revcomp = ''.join(comp[::-1])  
    
    StyI = 'CC[AT][AT]GG' 
    repl = 'NNNNNN' 
  
    pFwd = [f.start() for f in re.finditer(StyI, sequence)] 
    pRev = [r.start() for r in re.finditer(StyI, revcomp)]

    for i in pFwd:
        print(i, 'forward')

    for i in pRev:
        print(i, 'reverse')

    seqNew = re.sub(StyI,repl,sequence)
    print( "The sequence with StyI sites replaced is " + seqNew )

if args.e == 'styI':
  restriction_site_scan(sequence)







