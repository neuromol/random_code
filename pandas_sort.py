#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pandas Sort
■___■ 
(◕‿◕)
▐ __▐ 
.▆|▆.

"""
import pandas as pd 
import argparse
import sys

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Pandas sort using column names (compared to linux sort $columns)' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='Input file name ',required='True')
    parser.add_argument('-1', '--search1', help='Column name for first sort',required='True')
    parser.add_argument('-2', '--search2', help='Column name for second sort - leave out for None - default is None  ' ,default=None  )
    parser.add_argument('-o', '--output', help='outfile',required='True')
    results = parser.parse_args(args)
    return (results.input , results.search1, results.search2 , results.output)


def main(file, search1, search2, out ):
    
    input_file=pd.read_csv(file,sep="\t" , header=0)
    if search2 is None:
        sort_file = input_file.sort_values(by=[search1])    
    else:
        sort_file = input_file.sort_values(by=[search1,search2])    
    del input_file
    sort_file.to_csv(out,sep="\t",index=None)

if __name__ == '__main__':
    file, search1, search2, out  = check_arg(sys.argv[1:])

    main(file, search1, search2, out )
    
