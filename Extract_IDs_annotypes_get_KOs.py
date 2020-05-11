#This program accepts outputs of "get_gene_for_site_from_mlist.py" as inputs
#Uses annotations to fetch KEGG orthologies and pathway information for sites
#prints out global statistics
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
id_list=df1['Annotation'].fillna("").to_list()
atype_list=df1['Annotation'].fillna("").to_list()
#Parse accession IDs and infer annotation type
#Will be modified to use GFF file column #3 in the future
for i in range(0,len(id_list)):
    fields=id_list[i].split(';')
    if len(fields)>1:
        fields2=fields[0].split("_")
        if (len(fields2)==3):
            if (fields2[1]=="IGR"):
                atype_list[i]="IGR"
            elif (fields2[2]=="P"):
                atype_list[i]="promoter"
            elif (fields2[2]=="T"):
                atype_list[i]="terminator"
        elif (len(fields2)==2):
            atype_list[i]="gene"
        id_list[i]=fields[0].lstrip("ID=").rstrip("_P\n").rstrip("_T\n")
    elif (id_list[i]==""):
        atype_list[i] = "IGR"
id_series=pd.Series(id_list)
atype_series=pd.Series(atype_list)
#Add the IDs and Annotation Types to corresponding rows of dataframe
df2=pd.concat([df1,id_series,atype_series],axis=1).rename(columns={0:"ID",1:"Annotation Type"})
print (argstr.csvfile+"counts with duplicates...")
#Group sites by annotation type and report counts for each group
for name1,a_type in df2.groupby("Annotation Type"):
    print(name1 + ":" + str(a_type.count()))
print ("And without\n")
for name1,a_type in df2.drop_duplicates(subset=["chrom","start"]).groupby("Annotation Type"):
    print(name1 + ":" + str(a_type.count()))
ko_df=pd.read_csv(argstr.kos,sep=',',low_memory=False)
df3=pd.merge(df2,ko_df,on='ID',how='left')
if argstr.unique:
    df4=df3.sort_values(by=["ID","E-value"]).drop_duplicates(subset=["chrom","start","ID"])
    df4.to_csv(argstr.output,sep=',',index=False)
else:
    df3.to_csv(argstr.output, sep=',', index=False)
