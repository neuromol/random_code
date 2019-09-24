#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract rows from a list 
"""



import pandas as pd 
import argparse
import sys

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Merge two dataset together - eg extract a SNP list (Pandas merge and keep left index) ' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='Input file name ',required='True')
    parser.add_argument('-s', '--search', help='search file with list of elements to extract',required='True')
    parser.add_argument('-l', '--left', help='name of column on the search file  ' ,required='True'  )
    parser.add_argument('-r', '--right', help='name of the column on the large file',required='True' )
    parser.add_argument('-o', '--output', help='outfile',required='True')
    results = parser.parse_args(args)
    return (results.input , results.search, results.left , results.right , results.output)


def main(file,search, left,  right , out ):
    input_file=pd.read_csv(file,sep="\t" , header=0)
    search_file = pd.read_csv(search,sep="\t" , header=0)
    merged_file = search_file.merge(input_file, how="left" , left_on=left , right_on=right)
    del input_file
    del search_file
    merged_file.to_csv(out,sep="\t",index=None)


if __name__ == '__main__':
    file, search ,left,  right , out   = check_arg(sys.argv[1:])

    main(file, search, left,  right , out )