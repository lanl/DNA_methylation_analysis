#!/Users/shounak/anaconda3/bin/python3
#This program plots histograms to depict genome-wide methylation patterns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import matplotlib
import matplotlib.axes
matplotlib.rcParams['font.family']="monospace"
matplotlib.rcParams['font.monospace']="Courier New"
matplotlib.rcParams['font.size']=24
#argument handling
optparse = argparse.ArgumentParser()
optparse.add_argument("-c","--csvfile",help="list of methylation ratios")
optparse.add_argument("-t","--type",help="methlyation type:CpG, CHG or CHH")
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
aza_means=ratios.loc[:,(["Aza_1 "+argstr.date+" meth_ratio","Aza_2 "+argstr.date+" meth_ratio","Aza_3 "+argstr.date+" meth_ratio"])].loc[ratios['Type']==argstr.type].mean(axis=1).to_numpy()
#For control (untreated) samples
co_means=ratios.loc[:,(["Co_1 "+argstr.date+" meth_ratio","Co_2 "+argstr.date+" meth_ratio","Co_3 "+argstr.date+" meth_ratio"])].loc[ratios['Type']==argstr.type].mean(axis=1).to_numpy()

#means=pd.concat([aza_means,co_means],axis=1).rename(columns={0:"AZA",1:"CON"})
#hist_data=means.to_numpy()
#create a histogram for the 5-aza methylation calls...
plt.hist(aza_means,bins=np.arange(0.0,1,0.05),alpha=0.5,color="blue",label="AZA")
#... and the control methylation calls for a given methylation type (CpG,CHG or CHG)
plt.hist(co_means,bins=np.arange(0.0,1,0.05),alpha=0.5,color="red",label="CON")

#set the axis labels
plt.xlabel("Methylation ratio",fontsize=28)
plt.ylabel("Counts",fontsize=28)
#set the axis scales so we can compare plots
plt.xlim((0,1.0))
#plt.ylim((0,1.5E4))
#optional; tick label in scientific notation
plt.ticklabel_format(axis="y",scilimits=(2,4),useMathText=True)
# add the legend
plt.legend(fontsize=28,framealpha=1.0)
# and save the figure as a 300 DPI png file
plt.savefig(argstr.type+"_"+argstr.outfile+".png",dpi=300,format="png", bbox_inches='tight')
# close the plt object so that the above plots are not copied unintentionally,
# if this subroutine is called multiple times by the same parent python process
plt.close()

