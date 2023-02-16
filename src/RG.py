#! /usr/bin/env python
import argparse
import os
import sys
import gzip

def get_RGs(fastq_file):

    with gzip.open(fastq_file,'rb') as fin:        
        for line in fin:        
            line = line.strip().decode( "utf-8" )  
            break
        header = str(line).split(":")
        print(header)
        RG=header[0]+"_"+header[1]+"_"+header[2]+"_"+header[3]
        print(RG)
    SM = fastq_file.split("_")[0]
    LB = fastq_file.split("_")[1]
    LB = SM+"_"+LB 
    PU = RG+"."+LB 
    print("@RGID:",RG, "SM:", SM, "PL:Illumina", "LB:", LB, "PU:", PU)
    return RG, SM, LB, PU 

def run_picard(sam_file, outfile, RG, SM, LB, PU, PL): 
    I = "I=" + sam_file  
    O = "O=" + outfile 
    group = "RGID=" +RG+ " RGSM= " +SM +" RGPL=" + PL+ " RGLB="+ LB+ " RGPU=" +PU + " VALIDATION_STRINGENCY=SILENT" 
    cmd = "picard AddOrReplaceReadGroups " + I +" "+ O +" "+ group 
    os.system(cmd)     

def main(): 
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samfile',dest= "sam_file", 
                        required=True, help="sam file")
    parser.add_argument("-o", "--outfile", dest="outfile",required=True,
                        help="outfile")
    parser.add_argument("-f", "--fastq", dest="fastq_file", required=True,
                        help="fastq one pair file")
    parser.add_argument("-p", "--PL" , dest ="PL", default="Illumina") 
    args = parser.parse_args()
    fastq_file = args.fastq_file  
    sam_file = args.sam_file 
    outfile = args.outfile
    PL =args.PL 
    RG, SM, LB, PU = get_RGs(fastq_file)
    run_picard(sam_file, outfile, RG, SM, LB, PU, PL) 
if __name__ == "__main__":
    main()
