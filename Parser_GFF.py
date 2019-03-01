#!/usr/bin/python3
'''
Script that parses gff files and extracts position info for genes only
'''

import argparse

if  __name__== '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="input GFF file")
    args = parser.parse_args() 

file= args.f

#prepare file for data extraction
with open(file, 'r' ) as f:
   contents = f.read()

full_lines = contents.split('\n') #separate lines 

#create nested list where each item is a row of file
sep_lines=[]
for i in range(len(full_lines)):
   sep_lines.append(full_lines[i].split('\t'))

gff =[]
for item in sep_lines:
   if len(item) == 9:
      gff.append(item)

#parse gff file into a tuple containing chromosome #, gene start, gene end for each gene
gene_start=[]
gene_end = []
chrom_name= []
gene_id= [] #for gene name but idk what that is   

for i in range(len(gff)):
   row= gff[i] 
   start = int(row[3])  
   end = int(row[4])
   genex = row[2]
   chrom= row[0]
   id= row[8]
   #isolate genes only
   if genex == 'gene':
      chrom_name.append(chrom)
      gene_start.append(start)
      gene_end.append(end)
      gene_id.append(id)

            
list= [(chrom_name),(gene_start),(gene_end), (gene_id)]
gene_loc= zip(*list) #gff tuple of parsed information