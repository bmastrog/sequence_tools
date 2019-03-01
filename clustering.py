#!/usr/bin/python3

'''
Script that clusters differential expression data by hierarchical and k-means clustering
'''

import argparse
import pandas as p
from pandas import Series, DataFrame
from scipy import stats
from numpy import log2,mean
import scipy.cluster.hierarchy as h
from scipy.stats import zscore
from scipy.cluster.vq import kmeans,vq
import matplotlib as m
m.use('pdf')
from pylab import savefig

#set input file variable from command line
if  __name__== '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="input RNA seq data")
    args = parser.parse_args() 

file= args.f

# differential expression
DE = p.read_table(file ,delim_whitespace=True)
#DE = p.read_table("data2cluster.txt" ,delim_whitespace=True)
DEval = DE.values[:,1:7]  # root ctrl and treatment
DEval = log2(DEval.astype("double"))
p_val = Series([]) #capture t test
foldc = Series([]) #capture fold change
for i in range(DEval.shape[0]):
   logfc = mean(DEval[i,0:3]) - mean(DEval[i,3:6])
   tt = stats.ttest_ind(DEval[i,0:3],DEval[i,3:6])
   p_val[i] = tt[1]
   foldc[i] = abs(logfc) #abs(tt[0])

data_interest = DE[(p_val < 0.05) & (foldc > 3.0)]
genes = data_interest.values[:,0]
print (str(len(genes))+"\n")

# hierarchical clustering
data_whole = DE.values[:,1:]
data_whole = log2(data_whole.astype("double"))

x = h.linkage(data_whole, metric='correlation', method='average')
dendro = h.dendrogram(x)
savefig('dendrogram.pdf')
clstNum = h.fcluster(x, 4, criterion='maxclust')
print ("Sizes of hierachical clusters\n")
sum(clstNum == 1)
sum(clstNum == 2)
sum(clstNum == 3)
sum(clstNum == 4)

# kmeans clustering
from scipy.cluster.vq import kmeans2,vq
Z = zscore(data_whole, axis=1)
centroids,label = kmeans2(Z,4)
print ("Sizes of kmeans clusters\n")
sum(label == 0)
sum(label == 1)
sum(label == 2)
sum(label == 3)



