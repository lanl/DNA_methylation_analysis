#!/bin/python
#Take a CSV file with contig name and site coordinate
#and figure out if it lies within an annotated site
import argparse
import numpy as np
import pandas as pd
#Create object to parse arguments from the command line
optParse=argparse.ArgumentParser() #Create an argument parsing "object" called optParse
#define arguments and their flags
optParse.add_argument("-f","--csvfile",help="path to CSV input file")
optParse.add_argument("-m","--master",help="path to list of methylation sites")
optParse.add_argument("-o","--output",help="path to CSV output file")
#Parse arguments
argstr=optParse.parse_args()
#Read list of called methylation sites from CSV file; load into first dataframe
df1=pd.read_csv(argstr.csvfile,sep=",",low_memory=False)
#Load methylation site look up table into second dataframe
df2=pd.read_csv(argstr.master,sep=",",low_memory=False).rename(columns={"Scaffold":"chrom"}).drop_duplicates()
#Look up sites in first dataframe and fetch corresponding rows from second dataframe to merge; this is
#the validation; invalid sites are kept but receive no information from master lise (dataframe 2)
df3=pd.merge(df1,df2,on=["chrom","start"],how="left")
#write the merged data to the output file
df3.to_csv(argstr.output,index=False)
#delete the dataframes
del df1
del df2
del df3
#Inform the user of job completion
print ("Job complete")
