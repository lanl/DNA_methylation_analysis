#!/usr/bin/python
#pandas based program to build genewide methylation reports
#Shounak Banerjee 12-28-2019
#Steps are to load methylation reports with annotations into a dataframe
#Use built in 'group' method to split it into a list of dataframes,each corresponding to an annotation
#Then split by methylation type
#Fill all NaNs with 0's
#Gather reports
#Output format:
#Accession ID,Function,#CpGs Found,Mean CpG meth ratio,Max possible CpG,#CHGs Found,Mean CHG meth ratio,Max possible CHGs,#CHHs Found,Mean CHH meth ratio,
#Max Possible CHHs,Total C's methylated, Total C's methylable, Mean methylation (average of CpG, CHG and CHH meth ratios)


import numpy as np
import pandas as pd
import argparse
optParse=argparse.ArgumentParser() #Create an argument parsing "object" called optParse
optParse.add_argument("-f","--csvfile",help="path to CSV input file")
#optParse.add_argument("-g","--gff3file",help="path to GFF3 formatted annotation file")
optParse.add_argument("-o","--output",help="base name of CSV output files")
#optParse.add_argument("-i","--gff3out",help="path to GFF3 file with genes color coded by methylation value")
argstr=optParse.parse_args()
#We need to grab the name of column containing meth ratio data, typically with a timepoint prefix
report_df=pd.read_csv(argstr.csvfile,sep=',',low_memory=False).fillna(0.0) # create the dataframe from the CSV report
#report_df2=report_df.fillna(0.0) #replace all NaN's with 0's
meth_key=report_df.keys()[6]
report_df_grouped=report_df[["Annotation",meth_key]].groupby("Annotation") #First let's group the data by annotations
report_df_grouped_cpg=report_df.loc[report_df["Type"]=="CpG"][["Annotation",meth_key]].groupby("Annotation")
report_df_grouped_chg=report_df.loc[report_df["Type"]=="CHG"][["Annotation",meth_key]].groupby("Annotation")
report_df_grouped_chh=report_df.loc[report_df["Type"]=="CHH"][["Annotation",meth_key]].groupby("Annotation")
overall_report=["Accession,Mean "+meth_key+",Total "+meth_key+", Methylation sites available\n"]
cpg_report=["Accession,Type,Mean "+meth_key+",Total "+meth_key+", Methylation sites available\n"]
chg_report=["Accession,Type,Mean "+meth_key+",Total "+meth_key+", Methylation sites available\n"]
chh_report=["Accession,Type,Mean "+meth_key+",Total "+meth_key+", Methylation sites available\n"]
for i,j in report_df_grouped:
    df_temp=report_df_grouped.get_group(i)
    overall_report.append(str(i)+","+str(df_temp[meth_key].mean())+","+str(df_temp[meth_key].sum())+","+str(df_temp[meth_key].count())+"\n")
out_f = open(argstr.output, 'w')
for line in overall_report:
    out_f.write(line)
out_f.close()
del overall_report
for i, j in report_df_grouped_cpg:
    df_temp=report_df_grouped_cpg.get_group(i)
    cpg_report.append(str(i)+",CpG,"+str(df_temp[meth_key].mean())+","+str(df_temp[meth_key].sum())+","+str(df_temp[meth_key].count())+"\n")
out_cpg = open("CpG_"+argstr.output, 'w')
for line in cpg_report:
    out_cpg.write(line)
out_cpg.close()
del cpg_report
for i, j in report_df_grouped_chg:
    df_temp=report_df_grouped_chg.get_group(i)
    chg_report.append(str(i)+",CHG,"+str(df_temp[meth_key].mean())+","+str(df_temp[meth_key].sum())+","+str(df_temp[meth_key].count())+"\n")
out_chg = open("CHG_"+argstr.output, 'w')
for line in chg_report:
    out_chg.write(line)
out_chg.close()
del chg_report
for i, j in report_df_grouped_chh:
    df_temp=report_df_grouped_chh.get_group(i)
    chh_report.append(str(i)+",CHH,"+str(df_temp[meth_key].mean())+","+str(df_temp[meth_key].sum())+","+str(df_temp[meth_key].count())+"\n")
out_chh = open("CHH_"+argstr.output, 'w')
for line in chh_report:
    out_chh.write(line)
out_chh.close()
del chh_report
#gff_f=open(argstr.gff3file,'r')
#gff_out=open(argstr.gff3out,'w')
#out_f.write("Accession ID,Function,#CpGs Found,Mean CpG meth ratio,Max possible CpG,#CHGs Found,Mean CHG meth ratio,Max possible CHGs,#CHHs Found,Mean CHH meth ratio,Max Possible CHHs,Total C's methylated, Total C's methylable, Mean methylation (average of CpG, CHG and CHH meth ratios\n")
