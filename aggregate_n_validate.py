#!/Users/shounak/anaconda3/bin/python3
#This program supports the aggregator.py script
import numpy as np
import pandas as pd
import argparse
#argument handling
optparse = argparse.ArgumentParser()
optparse.add_argument("-c","--csvfile",help="list of methylation ratios")
optparse.add_argument("-t","--threshold",help="threshold mean values for both treated and untreated, to collect")
optparse.add_argument("-l","--lookup",help="look-up table to validate methylation sites")
optparse.add_argument("-d","--date",help="day in M/DD format, enclose within quotes")
optparse.add_argument("-o","--outfile",help="output histogram file basename")
argstr = optparse.parse_args()
#Read in the data
reads=pd.read_csv(argstr.csvfile,sep='\t',low_memory=False)
#Read in the validation table
val_tab=pd.read_csv(argstr.lookup,sep=',').rename(columns={"Scaffold":"chrom"}).drop_duplicates()
#Take the intersection of the reads and validation table to filter out the valid calls
ratios=pd.merge(reads,val_tab,on=["chrom","start"],how='inner')
#ratios.to_csv(argstr.outfile+"_all_validated.csv",sep=',',index=False)

#extract the relevant columns; need to replace this with generalized column list
#For 5-aza treated samples
aza_means=ratios.loc[:,(["Aza_1 "+argstr.date+" meth_ratio","Aza_2 "+argstr.date+" meth_ratio","Aza_3 "+argstr.date+" meth_ratio"])].mean(axis=1)
#For control (untreated) samples
co_means=ratios.loc[:,(["Co_1 "+argstr.date+" meth_ratio","Co_2 "+argstr.date+" meth_ratio","Co_3 "+argstr.date+" meth_ratio"])].mean(axis=1)

means=pd.concat([ratios,pd.concat([aza_means,co_means],axis=1).rename(columns={0:"AZA",1:"CON"})],axis=1).fillna("")
means.loc[(means["AZA"]>float(argstr.threshold)) & (means["CON"]>float(argstr.threshold))].to_csv(argstr.outfile,sep=',',index=False)