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
    parser = argparse.ArgumentParser(description='Make list of qQTLs rsids from gene list' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='input file',required='True')
    parser.add_argument('-k', '--key', help='list of genes', required='True')
    parser.add_argument('-o', '--output', help='output matrix', required='True')
    results = parser.parse_args(args)
    return (results.input , results.key, results.output )


def clean_me_plz(input_file, key_file , output):
    file = pd.read_csv(input_file , sep="\t", header=None)
    file.columns= ["rsid" , 'gene_id' , 't_stat']
    file_filt = file[(file["t_stat"] > 5 ) | (file["t_stat"] < -5 ) ]
    gene_list = pd.read_csv(key_file , sep="\t" , header=None)
    gene_list.columns = ["gene_id"]
    merged = gene_list.merge(file,on="gene_id")
    clean = merged["rsid"]
    clean.to_csv(output_file, sep="\t" , header=None , index=None)
if __name__ == '__main__':

    input_file, key_file , output_file = check_arg(sys.argv[1:])
    clean_me_plz(input_file, key_file , output_file)
