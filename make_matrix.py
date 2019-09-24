#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
▒▒▒▒▒▒▒▓
▒▒▒▒▒▒▒▓▓▓
▒▓▓▓▓▓▓░░░▓
▒▓░░░░▓░░░░▓
▓░░░░░░▓░▓░▓
▓░░░░░░▓░░░▓
▓░░▓░░░▓▓▓▓
▒▓░░░░▓▒▒▒▒▓
▒▒▓▓▓▓▒▒▒▒▒▓
▒▒▒▒▒▒▒▒▓▓▓▓
▒▒▒▒▒▓▓▓▒▒▒▒▓
▒▒▒▒▓▒▒▒▒▒▒▒▒▓
▒▒▒▓▒▒▒▒▒▒▒▒▒▓
▒▒▓▒▒▒▒▒▒▒▒▒▒▒▓
▒▓▒▓▒▒▒▒▒▒▒▒▒▓
▒▓▒▓▓▓▓▓▓▓▓▓▓
▒▓▒▒▒▒▒▒▒▓
▒▒▓▒▒▒▒▒▓

Make matrix 
"""


import pandas as pd
import sys
import argparse


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Make the matrix from files within a given directory' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='Input file',required='True' )
    parser.add_argument('-r', '--rsid', help='RSID list file to set the index on',required='True' )
    parser.add_argument('-t', '--tissue', help='tissue name file',required='True' )
    parser.add_argument('-o', '--output', help='output directory',required='True')
    parser.add_argument('-1', '--left', help='Name of column in RSID list',default="SNP")
    parser.add_argument('-2', '--right', help='Name of column in each tissue file',default="rs_id_dbSNP151_GRCh38p7")
    parser.add_argument('-n', '--name', help='name extention for each file eg ".clean for {tissue}.clean "',required='True')
    parser.add_argument('-s', '--stat', help='Name of the T-stat column or columns of interest',default="t_stat")
    results = parser.parse_args(args)
    return ( results.input  , results.rsid , results.tissue , results.output, results.left , results.right, results.name , results.stat)

    
       
def process_matrix(index_df, input_file , tissue , output, left , right, name , stat):
    ## where the magic happen loop each file to make a matrix
    for tiss in tissue:
        #print ("currently working on : " + tiss )
        file_label = input_file + tiss + name
        tissue_df = pd.read_csv(file_label, sep="\t", header=0)
        temp_df = tissue_df[[right , stat]].copy()
        new_column_name = tiss + "_T"
        temp_df.columns=[right , new_column_name]
        del tissue_df
        print ("printing the memory usage just for TS")
        print (temp_df.memory_usage(index=True).sum() )
        index_df = index_df.merge(temp_df, left_on=left , right_on=right, how="left")
        remove_rs_col = right + "_x"
        index_df = index_df.drop(right, axis=1)
      
        index_df = index_df.fillna("NA")
        del temp_df   
    index_df.to_csv(output, sep="\t" , index=None )
  
     

def get_tissues(tissue_file, left):
    ## returns list of tissues
    tissues_df = pd.read_csv(tissue_file,header=None)
    tissues_df.columns=['id']
    tissues_list = []
    for i in tissues_df.id :
        tissues_list.append(i)

    return tissues_list

def main(input_file  ,rsid , tissue , output, left , right, name , stat):
    tissues = get_tissues(tissue, left)
    index_df = pd.read_csv(rsid,sep="\t", header=0 )
    process_matrix(index_df, input_file , tissues , output, left , right, name , stat)


if __name__ == '__main__':
     input_file , rsid , tissue , output, left , right, name , stat = check_arg(sys.argv[1:])
     main(input_file , rsid , tissue , output, left , right, name , stat)

