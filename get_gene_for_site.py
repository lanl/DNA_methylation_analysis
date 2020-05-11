#!/bin/python
#Take a CSV file with contig name and site coordinate
#and figure out if it lies within an annotated site
import argparse
import numpy as np
import pandas as pd
def gff3_to_df(gff_file):
    gff_lines=gff_file.readlines()
    contig=[]
    start=[]
    annotation=[]
    for line in gff_lines[1:]:
        fields=line.split('\t')
        if fields[2]=="gene":
            for pos in range(int(fields[3]),int(fields[4])+1,1):
                contig.append(fields[0])
                start.append(pos)
                annotation.append(fields[8].rstrip('\n'))
    gff3_d={'Scaffold':contig,'start':start,'Annotation':annotation}
    del contig,start,annotation
    gff_df=pd.DataFrame(gff3_d)
    del gff3_d
    return gff_df
optParse=argparse.ArgumentParser() #Create an argument parsing "object" called optParse
optParse.add_argument("-f","--csvfile",help="path to CSV input file")
optParse.add_argument("-g","--gff3file",help="path to GFF3 file")
optParse.add_argument("-o","--output",help="path to CSV output file")
argstr=optParse.parse_args()
df1=pd.read_csv(argstr.csvfile) # open the input CSV file in read-only mode
with open(argstr.gff3file,'r') as gff3_f:
    df2=gff3_to_df(gff3_f)
gff3_f.close()
df3=pd.merge(df1,df2,on=["Scaffold","start"],how="left")
df3.to_csv(argstr.output,index=False)
del df1
del df2
del df3
print ("Job complete")
