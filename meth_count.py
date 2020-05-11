#!/usr/bin/python3
#Shounak Banerjee 10-27-2019
#Program to parse methylation calls in raw data CSVs, given methylation context
#Algorithm is to tally a call against a list of possible CpG, CHG and CHH methylation sites in a genome
#shortlists real methylation calls
#inputs: -f/--file <methylation call list in CSV format> -t/--type <CPG/CHH/CHG;case insensitive> -m/--methsites <methylation site list> -o/--output <output filename>
import argparse
import re
#Let's see what the user has given us

optParse=argparse.ArgumentParser() #Create an argument parsing "object" called optParse
optParse.add_argument("-f","--file",help="path to CSV input file")
optParse.add_argument("-t","--type",help="type of methylation; valid options are CpG,CHG and CHH")
optParse.add_argument("-m","--methsites",help="path to methylation site list")
optParse.add_argument("-o","--output",help="path to CSV output file")

argstr=optParse.parse_args() # store the arguments given/passed by the user

#Now, get to work and open the files (no error handling for now)

out_f=open(argstr.output,'w')# open the output CSV file in write-only mode

#Time to fetch the relevant list of methylation sites of a given type; no point in going through all the types!

meth_f=open(argstr.methsites,'r')
meth_f_lines=meth_f.readlines()
site_list=[]
for entry in meth_f_lines:
	if len(entry.split(","))>=4:
		if entry.split(",")[3]==argstr.type:
			site_list.append(entry)
meth_f.close()
csv_f=open(argstr.file,'r')
csv_list=csv_f.readlines()
out_f.write(meth_f_lines[0].replace("\n","")+","+csv_list[0].replace("\n","")+"\n")
csv_f.close()
#Now let's start reading the csv file, line by line
last_read=0
for line in site_list:
	#csv_f=open(argstr.file,'r') # open the input CSV file in read-only mode
	inp_fields=line.split(",")
	dataline=""
	for line2 in csv_list[last_read:]:
		csv_split=line2.split("\t")
		if inp_fields[2]=="+":
			if (inp_fields[0]==csv_split[0] and inp_fields[1]==csv_split[1]):
				dataline=str(csv_split[3])+","+str(csv_split[4])
				last_read=csv_list.index(line2)
				break
		if inp_fields[2]=="-":
			if (inp_fields[0]==csv_split[0] and inp_fields[1]==csv_split[2]):
				dataline=str(csv_split[3])+","+str(csv_split[4])
				last_read=csv_list.index(line2)
				break
	#csv_f.close()
	if line[len(line)-2]!=",":
		out_f.write(line.replace("\n","")+","+dataline.replace("\n","")+"\n")
	else:
		out_f.write(line.replace("\n","")+dataline.replace("\n","")+"\n")
out_f.close()
csv_f.close()
print("Methylation site list validation complete")

