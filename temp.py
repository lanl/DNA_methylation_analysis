# coding: utf-8
import numpy as np
import pandas as pd
import argparse
#Need to remove commas from KO Definitions so, open up the old file to process
old_ko=open("/Users/shounak/Documents/Algae_AOP/Epigenetics_CRS/Pico_KO_list_032820.csv",'r')
#and the new file to write to
new_ko=open("/Users/shounak/Documents/Algae_AOP/Epigenetics_CRS/Pico_KO_list_032820_regularized.csv",'w')
#Basically go line by line and only remove commas from the KO definitions if any are found
#They are replaced by spaces here but we can replace with any non-alphanumeric character
# that is NOT the FIELD SEPARATOR!
# Needed to do this for smooth reading into pandas dataframes

for line in old_ko.readlines():
    fields=line.split(",")
    new_ko.write(",".join(fields[0:5])+","+" ".join(fields[5:]).lstrip(" "))
#Now that we are done, close the file handles
new_ko.close()
old_ko.close()
#Load the reformatted KO table into a pandas DataFrame
df=pd.read_csv("/Users/shounak/Documents/Algae_AOP/Epigenetics_CRS/Pico_KO_list_032820_regularized.csv",sep=',')
#And load the mapping file
df2=pd.read_csv("/Users/shounak/Desktop/Mapped_KO_Pathway_list_w_descriptions.csv",sep=',')
#
import pandas as pd
import numpy as np
optParse=argparse.ArgumentParser() #Create an argument parsing "object" called optParse
optParse.add_argument("-f","--csvfile",help="path to CSV input file")
optParse.add_argument("-o","--output",help="path to CSV output file")
argstr=optParse.parse_args()
final_df=pd.read_csv("/Users/shounak/Documents/Algae_AOP/Epigenetics_CRS/Pico_KO_pathway_list.csv",sep=',')
hyper_df=pd.read_csv(argstr.csvfile,sep=',')
hyper_ko_pwy=pd.merge(hyper_df,final_df,on="ID",how='left')
#hyper_ko_pwy_unique=pd.merge(hyper_df,final_df,on="ID",how='left').drop_duplicates(subset="KO")
hyper_ko_pwy.to_csv(argstr.output+"all.csv",sep=',',index=False)
#hyper_ko_pwy.to_csv(argstr.output+"_unique.csv",sep=',',index=False)
