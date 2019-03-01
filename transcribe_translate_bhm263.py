#!/usr/bin/python3

"""
Script takes nucleotide sequences from input fasta file and runs seq tools functions 

"""
import argparse, seq_tools

if  __name__== '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="input fasta file w. nt sequences")
    args = parser.parse_args() 


sequence = seq_tools.load_FASTA(args.f) #set variable of input sequences 

seq_tools.get_reverse_complement(sequence)

seq_tools.get_transcript(sequence)

seq_tools.get_translation(sequence)




