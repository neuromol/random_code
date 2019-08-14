#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process Gtex qQTL
"""

import pandas as pd 
import numpy as np 
import argparse
import sys
import warnings
warnings.filterwarnings('ignore')

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Clean Gtex qQTL files' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='input file',required='True')
    parser.add_argument('-k', '--key', help='variant id to rsid from UK biobank', required='True')
    parser.add_argument('-g', '--gene', help='ensemble gene id to gene name', required='True')
    parser.add_argument('-o', '--output', help='output clean file', required='True')
    results = parser.parse_args(args)
    return (results.input , results.key, results.gene, results.output )


def main(input_file,key,gene,output):
    raw = pd.read_csv(input_file , sep="\t" , header=0)
    ## change the variant id in gtex so that we can link to rsid
    raw["variant_id"] =raw["variant_id"].str.replace("_", ":")
    ## get rid of the reference genome ad the end of the string
    raw["variant_id"] =raw["variant_id"].str.replace(":b37", "") 
    ## cal the Z-score
    raw["Z"] = raw["slope"] / raw["slope_se"]
    ## Edit the gene_id so that it gives the main gene ID (not isoform)
    raw["gene_id"] = raw["gene_id"].str.split(".").str[0]
    ##bring in the key file with has variant ID to rsid 
    key_file=pd.read_csv(key, sep="\t" , header=0)
    ## make a new dataframe with the info with want
    new = raw[["variant_id" , "gene_id" , "Z" ]]
    ### rename the columns so that variants is across both dataframes
    new.columns = ["variant" , "gene_id" , "Z"]
    ##bring in the gene map file
    gene_file=pd.read_csv(gene, sep="\t" , header=0)
    ## merge in to get the gene ID
    gene_combined = pd.merge(combined, gene_file, on="gene_id")
    ## The clean dataframe
    Final = gene_combined[["rsid" , "gene" , "Z"]]
    ###output 
    Final.to_csv(output, sep="\t" , index=None)


if __name__ == '__main__':

    input_file, key , gene, output = check_arg(sys.argv[1:])
    main(input_file, key , gene, output)
