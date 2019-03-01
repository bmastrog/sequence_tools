#!/usr/bin/python3
'''
Script that finds homologous pairs from blast outputs comparing two species
'''

import argparse

if  __name__== '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", help= "Species A:B .xml blast")
    parser.add_argument("-b", help= "Species B:A .xml blast")
    args = parser.parse_args() 

specAB = args.a #file comparing species a to b
specBA = args.b #file comparing species b to a


from Bio.Blast import NCBIXML
ab = NCBIXML.parse(open(specAB,"r"))
ba = NCBIXML.parse(open(specBA,"r"))
ab_list=[]
ba_list=[]

def hits(blast_file):
	"""This function would return high scoring pairs (hsp) from .xml blast output"""
	alist=[]
	for rec in blast_file:
		for alignment in rec.alignments:
			for hsp in alignment.hsps:
				if hsp.expect < 1e-20 : #sets a threshold value to e
					alist.append([rec.query,alignment.hit_def])
	
	final=[list(i) for i in set(map(tuple,alist))]
	return (final)

ab_list = hits(ab)
ba_list = hits(ba)

def reciprocal_hits(ab_list, ba_list):
	"""This function returns reciprocal hits ie putative orthologs from hsp of each xml output"""
	for x in ab_list: 
		if list([x[1],x[0]]) in ba_list:
			print (x[0] + "\t" + x[1])

reciprocal_hits(ab_list,ba_list)


