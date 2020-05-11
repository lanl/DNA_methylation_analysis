#!/usr/bin/python
import pandas as pd
df_25=pd.read_csv("7-25_all_validated_annotated.csv",sep=',',low_memory=False)
df_26=pd.read_csv("7-26_all_validated_annotated.csv",sep=',',low_memory=False)
result=pd.merge(df_25,df_26,how='inner',on=['chrom','start'])
df_26=result
df_27=pd.read_csv("7-27_all_validated_annotated.csv",sep=',',low_memory=False)
result=pd.merge(df_26,df_27,how='inner',on=['chrom','start'])
df_27=result
df_28=pd.read_csv("7-28_all_validated_annotated.csv",sep=',',low_memory=False)
result=pd.merge(df_27,df_28,how='inner',on=['chrom','start'])
df_28=result
df_31=pd.read_csv("7-31_all_validated_annotated.csv",sep=',',low_memory=False)
result=pd.merge(df_28,df_31,how='inner',on=['chrom','start'])
result.to_csv("Results_intersection.csv",sep=',',index=False)
result.loc[result['Type']=="CpG"].to_csv("Results_CpG_intersection.csv",sep=',',index=False)
result.loc[result['Type']=="CHG"].to_csv("Results_CHG_intersection.csv",sep=',',index=False)
result.loc[result['Type']=="CHH"].to_csv("Results_CHH_intersection.csv",sep=',',index=False)