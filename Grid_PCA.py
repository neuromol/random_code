#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PCA on grid
"""
import matplotlib
matplotlib.use('Agg')
import pandas as pd 
import argparse
import sys
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import numpy as np

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='PCA on the cluster (need the ram)' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='T stat file',required='True')
    parser.add_argument('-r', '--index', help='index column for the SNP/RS field', default="SNP" )
    parser.add_argument('-c', '--components', help='Number of components for PCA' ,required='True')
    parser.add_argument('-o', '--output', help='outfile of the PCs',required='True')
    results = parser.parse_args(args)
    return (results.input , results.index , results.components , results.output)

def main(file , index, comp , out):
    raw = pd.read_csv(file, sep="\t" , header=0, index_col=index)
    pca = PCA(n_components=comp)
    data = raw.replace("NA" , np.nan)
    data = raw.dropna()
    data = data.abs()
    data = StandardScaler().fit_transform(data)
    # print ("the length of dataframe is :" + str(len(data.columns)) ) 
    # print ("Printing datatypes; ")
    # print (data.dtypes)
    pca = PCA(n_components=int(comp))
    new_pca = pca.fit_transform(data)
    PCA_Stats = (pca.explained_variance_ratio_.cumsum() )
    # print ("Printing PCA explained variance stats per a PC ")
    print (PCA_Stats)
    new_pca = pd.DataFrame(new_pca , index=data.index)
    columns_list=[]
    for i in (range(0,len(new_pca.columns))):        
        dam_python = i + 1
        label = "PC" + str(dam_python)
        columns_list.append(label)
    new_pca.columns = columns_list
    new_pca.to_csv(out, sep="\t" )
    plt.scatter(new_pca.PC1 , new_pca.PC2 )
    plt.savefig("PCA.png" , dpi=600)

    

    
if __name__ == '__main__':
    file, index, comp , out  = check_arg(sys.argv[1:])

    main(file, index, comp,  out)
