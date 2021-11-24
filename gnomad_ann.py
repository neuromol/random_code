import numpy as np
import csv
import json
import sys
import pandas as pd
import requests
from pandas import json_normalize
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import os
import sys
import argparse
import time
from tqdm import tqdm

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Gnomad API for variants' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='Input file ',required='True')
    parser.add_argument('-o', '--output', help='output file', default='for_clump.txt')
    results = parser.parse_args(args)
    return (results.input , results.output) 

def get_gnomad_variants_for_gene(gene_id):
    query = """
    query Variants($gene_id: String!) {
      variant(variantId: $gene_id, dataset: gnomad_r3) {
        variantId
            genome{
      ac
      an
      homozygote_count
      af
    }
            in_silico_predictors {
      id
      value
    }
        sortedTranscriptConsequences {
        transcript_id
        consequence_terms
        canonical
        major_consequence
        polyphen_prediction
        sift_prediction
        lof
            }
          }
        }
        """
    variables = {
        "gene_id": gene_id,

    }
    
    headers = {"content-type": "application/json"}
    response = requests.post(
        "https://gnomad.broadinstitute.org/api", json={"query": query, "variables": variables}, headers=headers
    )

    response = json.loads(response.text)

    #errors = response.get("errors", [])
    # if errors:
    #     raise Exception(f"Error in response: {','.join(error['message'] for error in errors)}")

    variants = response
    return variants

def send_request1(gene_id):
#     gene_id = "ENSG00000010610"
    variants = get_gnomad_variants_for_gene(gene_id)
    df = json_normalize(variants)
    try: 
        df2 = pd.json_normalize(df['data.variant.in_silico_predictors'])
        temp = pd.DataFrame()
        for i in range(0,len(df2.columns)):
            temp = temp.append(pd.json_normalize(df2.iloc[:, i]))    
        moo = temp.T  
        new_header = moo.iloc[0]
        df_moo = moo[1:]
        df_moo.columns = new_header
        df_moo.reset_index(drop=True, inplace=True)
        df2 = pd.json_normalize(df['data.variant.sortedTranscriptConsequences'])
        df3 = pd.json_normalize(df2.iloc[:, 0])
        temp = pd.concat([df, df3, df_moo],1)
        temp = temp.drop(["data.variant.in_silico_predictors", "data.variant.sortedTranscriptConsequences"],1 ) 
        
        
    except:
        d = {'data.variant.variantId': gene_id, 'transcript_id': "." , 'consequence_terms' : "." , 'canonical':"." ,'major_consequence' : ".", 'polyphen_prediction': "." , 'sift_prediction':"." ,'lof': "." , "splice_ai": ".", 'revel': ".", 'primate_ai': "."  }
        df = pd.DataFrame(data=d, index=[0])
        return df
    return temp


def main(file_name , output_name): 
    input_file = pd.read_csv(file_name, sep="\t")
    input_file["CHROM"] = input_file['CHROM'].str[3:]
    input_file["gn_id"] = input_file.CHROM.astype(str) + "-" + input_file.POS.astype(str) + "-" + input_file.REF.astype(str) + "-" + input_file.ALT.astype(str)
    flat_list = input_file.gn_id.tolist()
    blank = pd.DataFrame()
    
    
    
    
    for snp in tqdm(flat_list):
        
        try :
            df = send_request1(snp)
            blank = blank.append(df)
            time.sleep(1)
        except:
            d = {'data.variant.variantId': snp, 'transcript_id': "." , 'consequence_terms' : "." , 'canonical':"." ,'major_consequence' : ".", 'polyphen_prediction': "." , 'sift_prediction':"." ,'lof': "." , "splice_ai": ".", 'revel': ".", 'primate_ai': "."  }
            df = pd.DataFrame(data=d, index=[0])
            blank = blank.append(df)
            time.sleep(1)


    merged = blank.merge(input_file, left_on="data.variant.variantId", right_on="gn_id")
    merged.to_csv(output_name, sep="\t", index=None)


if __name__ == '__main__':

    print ("""
    ┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘███┘┘┘┘┘┘┘
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘█████┘┘┘┘┘┘
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘███████████████┘┘┘┘┘
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘██████████████████┘┘┘┘
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘███████████████████┘██┘┘┘
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘██████████████┘┘███┘┘┘┘┘┘┘
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘██████████┘███┘┘┘█┘┘┘┘┘┘┘┘
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘██████┘┘┘┘┘███┘┘┘┘┘┘┘┘┘┘┘┘
┘┘┘┘┘┘┘┘┘┘┘┘██┘┘┘┘██████┘┘┘┘┘████┘┘┘┘┘┘┘┘┘┘┘
┘┘┘┘┘┘┘┘┘┘┘┘████┘┘┘██████┘┘┘┘┘████┘┘┘┘┘┘┘┘┘┘
┘┘┘┘┘┘┘┘┘┘┘┘┘█████┘███████┘┘┘┘┘████┘┘┘┘┘┘┘┘┘
┘┘┘┘┘┘┘┘┘┘┘┘███████┘██████┘┘┘┘┘┘┘███┘┘┘┘┘┘┘┘
┘┘┘┘┘┘┘┘┘┘┘┘┘█████████████┘┘┘┘┘┘┘┘████┘┘┘┘┘┘
┘┘┘┘██████┘┘┘┘████████████┘┘┘┘┘┘┘┘┘██┘┘┘┘┘┘┘
┘┘██████████┘┘┘███████████┘█████┘┘┘██┘┘┘┘┘┘┘
┘████┘┘┘┘┘███┘┘┘██████████████████████┘┘┘┘┘┘
███┘┘┘┘┘┘┘┘███┘┘┘█████████████████████┘┘┘┘┘┘
██┘┘┘█████┘┘███┘█████████████████████┘┘┘┘┘┘┘
██┘┘█████████████████████████████████┘┘┘┘┘┘┘
██┘┘█████████████████████████████████████┘┘┘
███┘┘┘██┘┘┘███┘┘██████████████████┘██┘┘████┘
████┘┘┘┘┘┘████┘┘┘████████████████┘┘██┘┘┘┘┘██
┘████████████┘┘┘██┘██████████┘┘┘┘┘┘██┘┘┘┘┘┘┘
┘┘┘┘███████┘┘┘┘┘┘┘┘┘┘████████┘┘┘┘┘┘██┘┘┘┘┘┘┘
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘███┘┘┘┘┘███████┘┘┘┘
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘█┘┘┘┘███████████┘┘
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘███┘┘██┘┘┘███┘
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘███┘┘┘██┘┘┘┘███
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘██┘┘┘┘████┘┘┘██
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘███┘┘██████┘┘┘██
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘██┘┘┘█████┘┘┘██
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘██┘┘┘┘███┘┘┘┘██
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘███┘┘┘┘┘┘┘┘┘██┘
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘████┘┘┘┘┘████┘
┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘┘█████████┘┘┘
___Gnomad api for a given variant___________
**** Please do not use this to flood the API - for intented use on a small variant list
    """)

    raw, output  = check_arg(sys.argv[1:])
    main(raw, output )