#!/usr/local/bin/python
import argparse
import re
optParse=argparse.ArgumentParser() #Create an argument parsing "object" called optParse
optParse.add_argument("-g","--gff3file",help="path to GFF3 formatted annotation file")
optParse.add_argument("-o","--output",help="path to GFF3 output file")
#Program to annotate 500 bp 5' & 3' UTRs as promoters and terminators
argstr=optParse.parse_args() # store the arguments given/passed by the user
gff3_lines=[]
p_start=0
p_end=0
p_name=""
p_contig=""
p_sense=''
t_start=0
t_end=0
t_name=""
t_contig=""
t_sense=''
with open (argstr.gff3file,'r') as gff3_f:
	gff3_lines=gff3_f.readlines()
gff3_contigs=[]
gff3_genes=[]
for line in gff3_lines[1:]:
	if (line.split('\t')[2]=="contig"):
		gff3_contigs.append(line)
	if (line.split('\t')[2]=="gene"):
		gff3_genes.append(line)
# First let's read the contig starts and ends into dictionaries
contig_starts={}
contig_ends={}
gff3_out=[]
gff3_out.append(gff3_lines[0])
for i in gff3_contigs:
	contig_starts[i.split('\t')[0]]=int(i.split('\t')[3])
	contig_ends[i.split('\t')[0]]=int(i.split('\t')[4])
#Handle first gene in gff3_genes
for line in gff3_genes:
	gff_line1=line.split('\t')
	if (gff_line1[6]=="+"):
		if (int(gff_line1[3])>501):
			p_start=int(gff_line1[3])-501
		else:
			p_start=1
		p_end=int(gff_line1[3])-1
		p_name=re.search("ID=.+?;",gff_line1[8]).group(0).replace(";","_P;")
		t_name=re.search("ID=.+?;",gff_line1[8]).group(0).replace(";","_T;")
		t_start=int(gff_line1[4])+1
		p_contig=gff_line1[0]
		t_contig=gff_line1[0]
		p_sense="+"
		t_sense="+"
		if (contig_ends[gff_line1[0]]-int(gff_line1[4])<500):
			t_end=contig_ends[gff_line1[0]]
		else:
			t_end=int(gff_line1[4])+500
	else:
		if (int(gff_line1[3])>501):
			t_start=int(gff_line1[3])-501
		else:
			t_start=1
		t_end=int(gff_line1[3])-1
		t_name=re.search("ID=.+?;",gff_line1[8]).group(0).replace(";","_P;")
		p_name=re.search("ID=.+?;",gff_line1[8]).group(0).replace(";","_T;")
		p_start=int(gff_line1[4])+1
		t_contig=gff_line1[0]
		p_contig=gff_line1[0]
		t_sense="-"
		p_sense="-"
		if (contig_ends[gff_line1[0]]-int(gff_line1[4])<500):
			p_end=contig_ends[gff_line1[0]]
		else:
			p_end=int(gff_line1[4])+500
	gff3_out.append(p_contig+"\tmanual\tgene\t"+str(p_start)+"\t"+str(p_end)+"\t.\t"+p_sense+"\t.\t"+p_name+"\n")
	gff3_out.append(t_contig+"\tmanual\tgene\t"+str(t_start)+"\t"+str(t_end)+"\t.\t"+t_sense+"\t.\t"+t_name+"\n")
with open(argstr.output,'w') as g_out:
	for line in gff3_out:
		g_out.write(line)