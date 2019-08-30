#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 12:22:43 2019

@author: josh
"""

import pandas as pd 
import numpy as np 
import argparse
import sys


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Make list of rsid from gene list' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='input file',required='True')
    parser.add_argument('-k', '--key', help='variant id to rsid from UK biobank', required='True')
    parser.add_argument('-g', '--gene', help='Gene list to extract', required='True')
    parser.add_argument('-o', '--output', help='output matrix', required='True')
    results = parser.parse_args(args)
    return (results.input , results.key, results.gene, results.output )



def clean_me_plz(input_file, key_file , gene_file, output):
    file = pd.read_csv(input_file , sep="\t", header=0)
    file["t_stat"] = file['slope'] / file['slope_se']
    file["variant_id"] = file["variant_id"].str[:-4]
    file["variant"] = file["variant_id"].str.replace('_',':')
    file['gene_id'] = file['gene_id'].str.split(".").str.get(0)
    map_file = pd.read_csv(key_file , sep="\t" , header= 0 )
    combined = file.merge(map_file,on="variant")
    file_filt = combined[(combined["t_stat"] > 5 ) | (combined["t_stat"] < -5 ) ]
    gene_list = pd.read_csv(gene_file, sep="\t" , header=None)
    gene_list.columns = ["gene_id"]
    merged = gene_list.merge(file_filt,on="gene_id")
    make_list = merged["rsid"].copy()
    make_list.to_csv(output, sep="\t" , header=None , index=None)


if __name__ == '__main__':

    input_file, key_file ,gene_file,  output_file = check_arg(sys.argv[1:])
    clean_me_plz(input_file, key_file ,gene_file, output_file)
