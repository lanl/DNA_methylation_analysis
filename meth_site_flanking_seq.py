#This python script reads a CSV formatted table of methylation sites
#and attaches, depending on the coordinate, 1.5 kb flanking regions
#numbers listed


import csv # module to read CSV files

import re # module to search for regular expressions in files; not in use now but will set up for sophisticated search and grabs

import argparse # module to handle or "parse" arguments

import numpy

def split_head_body(g_file,head_a,body_a):
	#a=[]
	#b=[]
	h_ct=-1
	for i in g_file:
		if (re.match('>',i)):
			h_ct=h_ct+1
			head_a.append(i.replace("\n","").replace(">",""))
			body_a.append("")
		else:
			body_a[h_ct]=body_a[h_ct]+i.replace("\n","")
	#head_a.append(a)
	#body_a.append(b)
	#print (headers)

#We will set up the following command line arguments
#-f --csvfile : path to CSV input file 
##-g --gff3file : path to GFF3 formatted annotation file
#-d --dnafile : path to FASTA formatted genome sequence file
#-o --output : path to CSV output file
#Let's see what the user has given us

optParse=argparse.ArgumentParser() #Create an argument parsing "object" called optParse
optParse.add_argument("-f","--csvfile",help="path to CSV input file")
optParse.add_argument("-s","--size",help="size of flanking region")
optParse.add_argument("-d","--dnafile",help="path to genome sequence file")
optParse.add_argument("-o","--output",help="path to CSV output file")

argstr=optParse.parse_args() # store the arguments given/passed by the user

#Now, get to work and open the files (no error handling for now)

csv_f=open(argstr.csvfile,'r') # open the input CSV file in read-only mode
out_f=open(argstr.output,'w')# open the output CSV file in write-only mode

#Now let's start reading the csv file, line by line

genome_f=open(argstr.dnafile,'r')
genome=genome_f.readlines()
head_a=[]
body_a=[]
split_head_body(genome,head_a,body_a)
flank_size=int(argstr.size)
#inp_csv_read=csv.reader(csv_f,dialect='excel')
#run a loop that iterates over all lines/rows in the CSV input
for line in csv_f:
	#store field values in an array
	inp_field=line.split(',')
	coord=int(inp_field[1])
	#Now, we know that the GenomeIDs are in Column 1
	#So we will use the first element of the array and search for matches in the
	#annotation file
	if (inp_field[0]!='' and inp_field!="Name"):
		if ((coord-flank_size)<1):
			coord_s=0
		else:
			coord_s=coord-flank_size-1
		if (coord+flank_size>len(body_a[head_a.index(inp_field[0])])):
			coord_e=len(body_a[head_a.index(inp_field[0])])
		else:
			coord_e=coord+flank_size
			seq=body_a[head_a.index(inp_field[0])][coord_s:coord_e]
			#print(genome_id.group(0).replace(";","")+","+seq+","+line)
			out_f.write(">"+inp_field[0]+":"+str(coord)+"\n"+seq+"\n\n")
csv_f.close()
out_f.close()
genome_f.close()
#count_fpkm.close()
quit()

	

