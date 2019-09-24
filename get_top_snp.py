#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Keep the top element - drop the rest
"""

import pandas as pd 
import argparse
import sys

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Drop duplicated elements (SNPs) based of a column on sorted files' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='Input file name ',required='True')
    parser.add_argument('-c', '--column', help='column name ',required='True')
    parser.add_argument('-o', '--output', help='outfile',required='True')
    results = parser.parse_args(args)
    return (results.input , results.column , results.output)


def main(file, column , out ):
    
    input_file=pd.read_csv(file,sep="\t" , header=0)
    clean = input_file.drop_duplicates([column]) 
    clean.to_csv(out,sep="\t",index=None)


if __name__ == '__main__':
    file, column , out  = check_arg(sys.argv[1:])
    main(file, column , out )
    