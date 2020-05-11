#!/usr/bin/python
#Shounak Banerjee 11-20-2019
#Program to aggregate methylation calls across a gene and report percent methylation
#inputs: -f/--file <methylation table> -g/--gff <gff3 file> -m/--methsites <methylation site list> -o/--output
#Output format:
#Accession ID,Function,#CpGs Found,Mean CpG meth ratio,Max possible CpG,#CHGs Found,Mean CHG meth ratio,Max possible CHGs,#CHHs Found,Mean CHH meth ratio,
#Max Possible CHHs,Total C's methylated, Total C's methylable, Mean methylation (average of CpG, CHG and CHH meth ratios)
import argparse
import matplotlib.pylab
from operator import itemgetter
def gather_gene_report(lines):
	cpg_tot=0
	cpg_found=0
	tot_cpg_meth_ratio=0.0
	mean_cpg_meth_ratio=0.0
	chg_tot=0
	chg_found=0
	tot_chg_meth_ratio=0.0
	mean_chg_meth_ratio=0.0
	chh_tot=0
	chh_found=0
	tot_chh_meth_ratio=0.0
	mean_chh_meth_ratio=0.0
	for site in lines:
		fields=site.split(",")
		if(fields[3]=="CpG"):
			cpg_tot=cpg_tot+1
			if fields[6]:
				if float(fields[6])>0.10:
					cpg_found=cpg_found+1
					tot_cpg_meth_ratio=tot_cpg_meth_ratio+float(fields[6])
		if cpg_tot>0:
			mean_cpg_meth_ratio=tot_cpg_meth_ratio/cpg_tot
		if(fields[3]=="CHG"):
			chg_tot=chg_tot+1
			if fields[6]:
				if float(fields[6])>0.10:
					chg_found=chg_found+1
					tot_chg_meth_ratio=tot_chg_meth_ratio+float(fields[6])
		if chg_tot>0:
			mean_chg_meth_ratio=tot_chg_meth_ratio/chg_tot
		if(fields[3]=="CHH"):
			chh_tot=chh_tot+1
			if fields[6]:
				if float(fields[6])>0.10:
					chh_found=chh_found+1
					tot_chh_meth_ratio=tot_chh_meth_ratio+float(fields[6])
		if chh_tot>0:
			mean_chh_meth_ratio=tot_chh_meth_ratio/chh_tot
	tot_c=cpg_tot+chg_tot+chh_tot
	meths_found=cpg_found+chg_found+chh_found
	mean_meth_ratio=mean_cpg_meth_ratio+mean_chg_meth_ratio+mean_chh_meth_ratio
	return (str(cpg_found)+","+str(mean_cpg_meth_ratio)+","+str(cpg_tot)+","+str(chg_found)+","+str(mean_chg_meth_ratio)+","+str(chg_tot)+","+str(chh_found)+","+str(mean_chh_meth_ratio)+","+str(chh_tot)+","+str(meths_found)+","+str(tot_c)+","+str(mean_meth_ratio))
#def write_color_coded_gff(line,meth_value):
#	red=int(meth_value*255)
#	green=128
#	blue=255-red
#	hex_code="#"+str("%02x%02x%02x" % (red,green,blue))
def write_color_coded_gff(line,meth_value,cmap)
	rgb = cmap(int(meth_value*255))[:3]
	return (line+";color="+plb.colors.tohex(rgb))
#start of main program
optParse=argparse.ArgumentParser() #Create an argument parsing "object" called optParse
optParse.add_argument("-f","--csvfile",help="path to CSV input file")
optParse.add_argument("-g","--gff3file",help="path to GFF3 formatted annotation file")
optParse.add_argument("-o","--output",help="path to CSV output file")
optParse.add_argument("-i","--gff3out",help="path to GFF3 file with genes color coded by methylation value")
argstr=optParse.parse_args()
try:
	csv_f=open(argstr.csvfile,'r') # open the input CSV file in read-only mode
except FileNotFoundError:
	print ("Error opening input CSV file:",argstr.csvfile) #if something goes wrong; inform the user and
	quit() #quit the program
gff_f=open(argstr.gff3file,'r')
gff_lines=gff_f.readlines()
gff_f.close()
gff_genes=[]
gff_out=open(argstr.gff3out,'w')
gff_out.write("##gff-version 3\n")
for line in gff_lines[1:]:
	if (len(line.split('\t'))==9):
		if (line.split('\t')[2]=="gene"):
			gff_genes.append(tuple(line.split('\t')))
		else:
			gff_out.write(line)
sorted_gff_gene_tuples=sorted(gff_genes,key=itemgetter(0,3))
sorted_gff_genes=[]
for i in sorted_gff_gene_tuples:
	sorted_gff_genes.append("\t".join(i))
del sorted_gff_gene_tuples
del gff_genes
acc_list=[]
out_f=open(argstr.output,'w')
out_f.write("Accession ID,Function,#CpGs Found,Mean CpG meth ratio,Max possible CpG,#CHGs Found,Mean CHG meth ratio,Max possible CHGs,#CHHs Found,Mean CHH meth ratio,Max Possible CHHs,Total C's methylated, Total C's methylable, Mean methylation (average of CpG, CHG and CHH meth ratios\n")
csv_lines_pre=csv_f.readlines()[1:]
csv_f.close()
csv_lines_pre2=[]
for a_line in csv_lines_pre:
	csv_lines_pre2.append(tuple(a_line.split(',')))
del csv_lines_pre
csv_line_tuples=sorted(csv_lines_pre2,key=itemgetter(0,1))
csv_lines=[]
for i in csv_line_tuples:
	csv_lines.append(",".join(i))
del csv_line_tuples
#phew, we just sorted those all the data!
last_read=0
my_cmap=plb.get_cmap("magma",1)
for line2 in sorted_gff_genes:
	gff_out.write(line2)
	accession=line2.split("\t")[8].split(";")[0].lstrip("ID=")+";"
	feat=line2.split("\t")[2]
	acc_list.append(accession)
	relevant_lines=[]
	loop_breaker=""
	func=""
	for line3 in csv_lines[last_read:]:
		last_read=last_read+1
		if (accession in line3):
			loop_breaker=accession
			relevant_lines.append(line3)
			func=line3.split(',')[4]
		if loop_breaker!="" and loop_breaker not in line3:
			break
	report=gather_gene_report(relevant_lines)
	meth_v=float(report.split(",")[11])
	out_f.write(accession+","+func+","+report+"\n")
	gff_out.write(write_color_coded_gff(line2.replace("\t"+feat+"\t","\tmethylated_C\t"),meth_v,my_cmap)+"\n")
gff_out.close()
out_f.close()
