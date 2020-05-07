import myvariant
import pandas as pd
import sys
import argparse
mv = myvariant.MyVariantInfo()

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Get nearby variants filtering by CADD and size ' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--cadd', help='above a CADD score',required='True')
    parser.add_argument('-i', '--input', help='Input SNP list with header',required='True')
    parser.add_argument('-s', '--size', help='size in bases',required='True')
    parser.add_argument('-o', '--output', help='output file',required='True')
    results = parser.parse_args(args)
    return (results.cadd, results.input, results.size, results.output)


def find_miss(val_pos,chrom,size, cadd):
    ''' val_pos the position of variant (hg19) , chrom = chromosome 15 add the chr1 , size = window size'''
    ### get the first bit of the commands out of the way
    #snpeff_command = " AND snpeff.ann.effect:" + snpeff
    cadd = " AND cadd.phred:>" + str(cadd)
    #gnomad_freq = " AND gnomad_genome.af.af:>" + str(freq)
    ### set the sizes
    start_pt = val_pos - size
    if start_pt < 0 :
        start_pt = 0 
    chrom = "chr" + str(chrom)
    finish_pt = val_pos + size 
    ##### make the command AND _exists_:dbsnp
    start_command = chrom + ":" + str(start_pt) + "-" + str(finish_pt) 
    #query =  start_command + cadd + snpeff_command + gnomad_freq
    query =  start_command + cadd
    testing = mv.query(query, fields='snpeff.ann.effect, cadd.phred, snpeff.ann.genename, dbsnp.rsid, cadd.polyphen.cat, cadd.sift.cat, gnomad_genome.af.af, clinvar.gene.id, vcf , gnomad_exome.af.af' , as_dataframe=True, size=1000, assembly='hg19')
    print (len(testing))
    try:
        First = testing['snpeff.ann'].apply(pd.Series)[0]
        First = First.apply(pd.Series)
        moocow = testing.drop(["_id", "_score", "cadd._license", "snpeff._license","gnomad_exome._license", "snpeff.ann","gnomad_genome._license","dbsnp._license"],axis=1)
        
        moocow = pd.concat([moocow, First], axis=1)
        return moocow
    except:
        return 
    
def process_input(input_file):
    demo = pd.DataFrame()
    for chunk in pd.read_csv(input_file, header=0, chunksize=1000):
        chunk.columns = ["SNP"]
        chunk_list = chunk.SNP.to_list()
        temp= mv.querymany(chunk_list, scopes='dbsnp.rsid', fields='vcf, chrom, dbsnp.rsid',  as_dataframe=True, size=1000, assembly='hg19')
        demo = demo.append(temp)
    demo = demo.drop_duplicates("dbsnp.rsid")
    return demo

def main(cadd , input_file , size , output_file):   
    input_df = process_input(input_file)
    moo = pd.DataFrame() 
    for index, row in input_df.iterrows():
        chrom = row['chrom']
        pos = int(row['vcf.position']) 
        temp = find_miss(pos, chrom , int(size), cadd)
        moo = moo.append(temp)
    moo.to_csv(output_file,sep="\t", index=None)
    
if __name__ == '__main__':
    print ("""
    Get variants near centroids
    -----------------------------------
    
    
    """)
    
    cadd , input_file , size , output_file = check_arg(sys.argv[1:])
    
    print ("cadd score : " + str(cadd))
    print ("input file : " + input_file)
    print ("output file : " + output_file)
    print ("window size: " + str(size) + " bases")
    
    main(cadd , input_file , size , output_file )    
    
