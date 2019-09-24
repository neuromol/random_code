#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:14:43 2019

@author: josh
"""

import matplotlib
matplotlib.use('Agg')
import pandas as pd 
import argparse
import sys

import matplotlib.pyplot as plt
from sklearn.decomposition import IncrementalPCA
import numpy as np


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='IncrementalPCA on the cluster (need the ram)' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='T stat file',required='True')
    parser.add_argument('-c', '--components', help='Number of components for PCA' )
    parser.add_argument('-o', '--output', help='outfile of the PCs',required='True')
    results = parser.parse_args(args)
    return (results.input , results.components , results.output)

def main(file , comp , out):
    raw = pd.read_csv(file, sep="\t" , header=0, index_col="rs_id_dbSNP151_GRCh38p7" , chunksize=100000)
    index_list=[]
    ipca = IncrementalPCA(n_components=int(comp))
    for chunk in raw:       
        ipca.partial_fit(chunk)
        for b in chunk.index.tolist() :
            index_list.append(b)
    mean = ipca.mean_
    stddev = np.sqrt(ipca.var_)
    print (mean)
    print (stddev)
    Xtransformed = None
    second_loop = pd.read_csv(file, sep="\t" , header=0, index_col="rs_id_dbSNP151_GRCh38p7" , chunksize=100000)    
    output_df = pd.DataFrame()
    for chunk in second_loop:
        moo = ipca.transform(chunk)       
        temp = pd.DataFrame(moo)
        output_df = output_df.append(temp)
            
    columns_list=[]
    for i in (range(0,len(output_df.columns))):        
        dam_python = i + 1
        label = "PC" + str(dam_python)
        columns_list.append(label)
    
    print (columns_list)
    output_df.columns= columns_list
    output_df.index = index_list
    print (output_df)
    #print (index_list)
    PCA_Stats = (ipca.explained_variance_ratio_.cumsum() )
    PCA_stats_df = pd.DataFrame(PCA_Stats)
    print (PCA_stats_df)
    output_df.to_csv(out, sep="\t" )
    plt.scatter(output_df.PC1 , output_df.PC2 )
    plt.savefig("IPCA.png" , dpi=600)
    
if __name__ == '__main__':
    file, comp , out  = check_arg(sys.argv[1:])

    main(file, comp,  out)
