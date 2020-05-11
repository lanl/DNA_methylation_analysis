#!/bin/python
#Shounak Banerjee 10-27-19
#Program to enumerate all possible methylation sites in a genome
#Need to find CpG CG, CHG C[ATC]G and CHH C[ATC][ATC]
#
#
import re
import argparse
def find_CpG(name,seq):#outputs a string of the fasta section header and start coordinate of a CpG methylation site 
	#fortunately for CpG's, the forward and reverse sites are 'adjacent'
	cpg_site=re.finditer("[Cc][Gg]",seq)
	output=""
	if cpg_site:
		for site in cpg_site:
			output=output+name+","+str(site.span()[0]+1)+",+,CpG\n"+name+","+str(site.span()[1])+",-,CpG\n"
	return output
def find_CHG(name,seq):#outputs a string of the fasta section header and start coordinate of a CHG methylation site 
	#first let's handle the forward strand
	chg_site=re.finditer("[Cc][AaTtCc][Gg]",seq)
	output=""
	if chg_site:
		for site in chg_site:
			output=output+name+","+str(site.span()[0]+1)+",+,CHG\n"
	#and because we can't forget the reverse strand
	chg_site_r=re.finditer("[Cc][AaTtGg][Gg]",seq)
	if chg_site_r:
		for site in chg_site_r:
			output=output+name+","+str(site.span()[1])+",-,CHG\n"
	return output
def find_CHH(name,seq):#outputs a string of the fasta section header and start coordinate of a CHH methylation site
	#first let's handle the forward strand
	chh_site=re.finditer("[Cc][AaTtCc][AaTtCc]",seq)
	output=""
	if chh_site:
		for site in chh_site:
			output=output+name+","+str(site.span()[0]+1)+",+,CHH\n"
	chh_site_r=re.finditer("[AaTtGg][AaTtGg][Gg]",seq)
	if chh_site_r:
		for site in chh_site_r:
			output=output+name+","+str(site.span()[1])+",-,CHH\n"
	return output
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
#Main program
optParse=argparse.ArgumentParser() #Create an argument parsing "object" called optParse
optParse.add_argument("-d","--dnafile",help="path to genome sequence file")
optParse.add_argument("-o","--output",help="path to CSV output file")

argstr=optParse.parse_args() # store the arguments given/passed by the user

#Now, get to work and open the files (no error handling for now)

out_f=open(argstr.output,'w')# open the output CSV file in write-only mode

#Now let's start reading the csv file, line by line

genome_f=open(argstr.dnafile,'r')
genome=genome_f.readlines()
head_a=[]
body_a=[]
split_head_body(genome,head_a,body_a)
output_collector=""
out_f.write("Scaffold,start,Strand,Type\n")
for i in range(0,len(head_a)):
	out_f.write(find_CpG(head_a[i],body_a[i])+find_CHG(head_a[i],body_a[i])+find_CHH(head_a[i],body_a[i]))
out_f.close()
