#!/usr/bin/python3

"""
This script contains functions for DNA sequence analysis -- parses sequences from FASTA files;
provides reverse complement, mRNA transcript, and translated protein sequence.
"""

def load_fasta(filename):
   "This function isolates the nucleotide sequences from a FASTA file and returns them"
   with open(filename,'r') as f:
      seq = []
      for line in f:    
         if line.startswith('>'): #Skips annotation (start with >)
            pass
         else:
            add = line.rstrip('\n') #removes newline symbols
            seq.append(add) #adds the isolated sequence(s) to a list     
   return seq
 
 
def load_FASTA(filename):
   with open(filename) as file:
      entries = file.read().split('>')[1:] #read the file and separate entries that start with '>'
      partitioned_entries = [entry.partition('\n') for entry in entries] # make 3-item tuple with header, newline, sequence
      pairs = [(entry[0], entry[2]) for entry in partitioned_entries]    # separate header line and sequence
      sequence = [(pair[1].replace('\n', '')) for pair in pairs]  # remove newlines in sequence & save only sequence
   return sequence  
  
valid= "ACTGN" #sets the allowed nucleotides for sequence
   
def get_reverse_complement(list):
   "This function returns the reverse complement of a nucleotide sequence"
   for item in list:
      if (all(i in valid for i in item) == True): #checks sequence for allowed nt only
         intab = "ACTGN"  
         outtab = "TGACN"
         trantab = item.maketrans(intab, outtab)
         print("The reverse complement is:" , item.translate(trantab) ) #adds translated seqs to a list
      else:
         print ("Sequence does not contain only nucleotides") 
        
   
    


def get_transcript(list):
   "This function returns the mRNA transcript from DNA"
   for item in list:
      if (all(i in valid for i in item) == True): #checks sequence for allowed nt only
         intab = "ACTGN"
         outtab = "ACUGN"
         trantab = item.maketrans(intab, outtab)
         print("The mRNA transcript is:" , item.translate(trantab) )
      else:
         print ("Sequence does not contain only nucleotides") 
       
   


def get_translation(list):
    "This translates a DNA sequence to a protien sequence"
    #Set dictionary of genetic code for translate function
    genetic_code = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
    'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
    'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
    'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
    'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
    'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
    'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
    'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
    'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
    'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
    'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
    'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGA': 'R',
    'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
    'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
    'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
    'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}
    protein=''

    for item in list:
       if (all(i in valid for i in item) == True):
          for n in range(0,len(item),3): #for every 3 bases over length of string
              if item[n:n+3] in genetic_code: #if letters make up a codon, return the protein
                  codon=item[n:n+3]
                  protein += genetic_code[ codon ]
          print( "The translated protien sequence is:", protein)
       else:
         print ("Sequence does not contain only nucleotides") 



#test sequences 
if __name__ == '__main__':
   get_transcript(["ATATDSTAV"])
   get_transcript(["ATATTGTAC"])
   get_reverse_complement(["ATATTGTAC"])
   get_reverse_complement(["ATATDSTAV"])
   get_translation(["ATATTGTAC"])
   get_translation(["ATATDSTAV"])




