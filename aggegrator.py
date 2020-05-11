#!/usr/bin/python
import numpy as np
import pandas as pd
df_25=pd.read_csv("7-25/CHG_result_table.txt_agt20.txt",sep='\t',low_memory=False)
df_26=pd.read_csv("7-26/CHG_result_table.txt_agt20.txt",sep='\t',low_memory=False)
#Retain only sites that appear on 7-25 AND 7-26
result=pd.merge(df_25,df_26,how='inner',on=['chrom','start'])
df_26=result
df_27=pd.read_csv("7-27/CHG_result_table.txt_agt20.txt",sep='\t',low_memory=False)
#Retain only sites that appear on 7-25 AND 7-26 AND 7-27
result=pd.merge(df_26,df_27,how='inner',on=['chrom','start'])
df_27=result
df_28=pd.read_csv("7-28/CHG_result_table.txt_agt20.txt",sep='\t',low_memory=False)
#Retain only sites that appear on 7-25 AND 7-26 AND 7-27 AND 7-28
result=pd.merge(df_27,df_28,how='inner',on=['chrom','start'])
df_28=result
df_31=pd.read_csv("7-31/CHG_result_table.txt_agt20.txt",sep='\t',low_memory=False)
#Retain only sites that appear on 7-25 AND 7-26 AND 7-27 AND 7-28 AND 7-31
result=pd.merge(df_28,df_31,how='inner',on=['chrom','start'])
#Write the resulting time series to a CSV file
result.to_csv("CHG_results_intersection.tsv",sep='\t')
#Repeat the process for CHH calls
df_25=pd.read_csv("7-25/CHH_result_table.txt_agt20.txt",sep='\t',low_memory=False)
df_26=pd.read_csv("7-26/CHH_result_table.txt_agt20.txt",sep='\t',low_memory=False)
result=pd.merge(df_25,df_26,how='inner',on=['chrom','start'])
df_26=result
df_27=pd.read_csv("7-27/CHH_result_table.txt_agt20.txt",sep='\t',low_memory=False)
result=pd.merge(df_26,df_27,how='inner',on=['chrom','start'])
df_27=result
df_28=pd.read_csv("7-28/CHH_result_table.txt_agt20.txt",sep='\t',low_memory=False)
result=pd.merge(df_27,df_28,how='inner',on=['chrom','start'])
df_28=result
df_31=pd.read_csv("7-31/CHH_result_table.txt_agt20.txt",sep='\t',low_memory=False)
result=pd.merge(df_28,df_31,how='inner',on=['chrom','start'])
result.to_csv("CHH_results_intersection.tsv",sep='\t')
#...and CpG calls
df_25=pd.read_csv("7-25/CpG_result_table.txt_agt20.txt",sep='\t',low_memory=False)
df_26=pd.read_csv("7-26/CpG_result_table.txt_agt20.txt",sep='\t',low_memory=False)
result=pd.merge(df_25,df_26,how='inner',on=['chrom','start'])
df_26=result
df_27=pd.read_csv("7-27/CpG_result_table.txt_agt20.txt",sep='\t',low_memory=False)
result=pd.merge(df_26,df_27,how='inner',on=['chrom','start'])
df_27=result
df_28=pd.read_csv("7-28/CpG_result_table.txt_agt20.txt",sep='\t',low_memory=False)
result=pd.merge(df_27,df_28,how='inner',on=['chrom','start'])
df_28=result
df_31=pd.read_csv("7-31/CpG_result_table.txt_agt20.txt",sep='\t',low_memory=False)
result=pd.merge(df_28,df_31,how='inner',on=['chrom','start'])
result.to_csv("CpG_results_intersection.tsv",sep='\t')