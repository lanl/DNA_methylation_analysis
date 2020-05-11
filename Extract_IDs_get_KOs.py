#This program accepts outputs of "get_gene_for_site_from_mlist.py" as inputs
#Uses annotations to fetch KEGG orthologies and pathway information for sites
import numpy as np
import pandas as pd
import argparse
optParse=argparse.ArgumentParser() #Create an argument parsing "object" called optParse
optParse.add_argument("-f","--csvfile",help="path to CSV input file")
optParse.add_argument("-k","--kos",help="path to list of KOs")
optParse.add_argument("-o","--output",help="path to CSV output file")
optParse.add_argument("-u","--unique",help="If True, then keep only highest confidence KO/pathway assignment for a row in input",default=False)
argstr=optParse.parse_args()
df1=pd.read_csv(argstr.csvfile,sep=",",low_memory=False)
#Extract the Annotation column in the methylation site list
id_list=df1['Annotation'].fillna("").to_list()
for i in range(0,len(id_list)):
    fields=id_list[i].split(';')
    if len(fields)>1:
        #Remove _P and _T extensions to accession IDs to enable look up in accession-->KO map
        id_list[i]=fields[0].lstrip("ID=").rstrip("_P\n").rstrip("_T\n")
id_series=pd.Series(id_list)
#Add the ID's as a separate column to df1 to create df2
df2=pd.concat([df1,id_series],axis=1).rename(columns={0:"ID"})
ko_df=pd.read_csv(argstr.kos,sep=',',low_memory=False)
#Assign the KOs and Pathways
df3=pd.merge(df2,ko_df,on='ID',how='left')
#If user chooses to turn on the 'unique' option
if argstr.unique:
    df4=df3.sort_values(by=["ID","E-value"]).drop_duplicates(subset=["chrom","start","ID"])
    df4.to_csv(argstr.output,sep=',',index=False)
else:
    df3.to_csv(argstr.output, sep=',', index=False)
