#!/Users/shounak/anaconda3/bin/python3
import argparse
import re
from operator import itemgetter
from Bio import SeqIO
from Bio import SeqUtils
from Bio import Seq

# To label the IGRs we will first create annotation contigs
# Then we will take the genomic scaffold/contigs and associated annotations
# Calculate the "complement set" ie. stuff not annotated
# Need to have a global counter to assign IGR accession IDs

optParse = argparse.ArgumentParser()  # Create an argument parsing "object" called optParse
optParse.add_argument("-f","--fasta",help="path to genome sequence file")
optParse.add_argument("-g", "--gff3file", help="path to GFF3 formatted annotation file")
optParse.add_argument("-o", "--output", help="path to GFF3 output file")
# Program to annotate 500 bp 5' & 3' UTRs as promoters and terminators
argstr = optParse.parse_args()  # store the arguments given/passed by the user
gff3_lines = list()
p_start = 0
p_end = 0
p_name = ""
p_contig = ""
p_sense = ''
t_start = 0
t_end = 0
t_name = ""
t_contig = ""
t_sense = ''
with open(argstr.gff3file, 'r') as gff3_f:
	gff3_lines = gff3_f.readlines()
#gff3_contigs = []
gff3_genes = []
for line in gff3_lines[1:]:
	if (line != "\n" and len(line.split('\t'))==9):
		#if (line.split('\t')[2] == "contig"): #03-24-20 need to refactor this to rely on genome sequence file
		#	gff3_contigs.append(line)
		if (line.split('\t')[2] == "gene"):
			gff3_genes.append(line)
# First let's read the contig starts and ends into dictionaries
contigs = {}
gff3_out = list()
gff3_out.append(gff3_lines[0])
genome_seqs=SeqIO.parse(argstr.fasta,"fasta")
#for i in gff3_contigs:
for i in genome_seqs:
	fields1 = i.id
	contigs[fields1] = (1, len(i.seq))
	#gff3_out.append(i)
print (contigs)
# Handle first gene in gff3_genes
for line in gff3_genes:
	gff3_out.append(line)
	gff_line1 = line.split('\t')
	if (gff_line1[6] == "+"):
		if (int(gff_line1[3]) > 501):
			p_start = int(gff_line1[3]) - 501
		else:
			p_start = 1
		p_end = int(gff_line1[3]) - 1
		p_name = re.search("ID=.+?;", gff_line1[8]).group(0).replace(";", "_P;")
		t_name = re.search("ID=.+?;", gff_line1[8]).group(0).replace(";", "_T;")
		t_start = int(gff_line1[4]) + 1
		p_contig = gff_line1[0]
		t_contig = gff_line1[0]
		p_sense = "+"
		t_sense = "+"
		if (contigs[gff_line1[0]][1] - int(gff_line1[4]) < 500):
			t_end = contigs[gff_line1[0]][1]
		else:
			t_end = int(gff_line1[4]) + 500
	else:
		if (int(gff_line1[3]) > 501):
			t_start = int(gff_line1[3]) - 501
		else:
			t_start = 1
		t_end = int(gff_line1[3]) - 1
		t_name = re.search("ID=.+?;", gff_line1[8]).group(0).replace(";", "_P;")
		p_name = re.search("ID=.+?;", gff_line1[8]).group(0).replace(";", "_T;")
		p_start = int(gff_line1[4]) + 1
		t_contig = gff_line1[0]
		p_contig = gff_line1[0]
		t_sense = "-"
		p_sense = "-"
		if (contigs[gff_line1[0]][1] - int(gff_line1[4]) < 500):
			p_end = contigs[gff_line1[0]][1]
		else:
			p_end = int(gff_line1[4]) + 500
	gff3_out.append(
		p_contig + "\tmanual\tgene\t" + str(p_start) + "\t" + str(p_end) + "\t.\t" + p_sense + "\t.\t" + p_name + "\n")
	gff3_out.append(
		t_contig + "\tmanual\tgene\t" + str(t_start) + "\t" + str(t_end) + "\t.\t" + t_sense + "\t.\t" + t_name + "\n")
anno_list = []
for item in gff3_out[1:]:
	fields2 = item.split('\t')
	if (len(fields2) == 9):
		if (fields2[2] == "gene"):
			anno_list.append((fields2[0], int(fields2[3]), int(fields2[4])))
sorted_anno_list = sorted(anno_list, key=itemgetter(0, 1))  # then sort by contig name
contig1 = ''
start = 0
end = 0
contig2 = ''
start2 = 0
end2 = 0
igr_count = 0
igr_base = 'NSC_IGR_'
igr_tuple = []
index_anno = 0
anno_len = len(sorted_anno_list) - 1
index_anno2 = 1
while (index_anno < anno_len-1):
	contig1 = sorted_anno_list[index_anno][0]
	start = sorted_anno_list[index_anno][1]
	end = sorted_anno_list[index_anno][2]
	index_anno2=index_anno+1
	if (index_anno == 0 and start > 1):
		igr_count += 1
		igr_tuple.append((contig1, (igr_base + str(igr_count)), 1, start - 1))
	print("index 1:"+str(index_anno)+"\tindex 2:"+str(index_anno2)+"\n")
	while(index_anno2 < anno_len):
		contig2 = sorted_anno_list[index_anno2][0]
		start2 = sorted_anno_list[index_anno2][1]
		end2 = sorted_anno_list[index_anno2][2]
		if contig1 == contig2 and index_anno2!=anno_len-1:
			if (start2 <= end + 1):
				end = end2
				index_anno2 += 1
			else:
				index_anno += index_anno2+ 1
				#print("Non-overlap found:" + str(end) + ":" + str(start2))
				igr_count += 1
				igr_tuple.append((contig1, (igr_base + str(igr_count)), end + 1, start2 - 1))
				break
		elif contig1!=contig2 and index_anno2!=anno_len-1:
			index_anno = index_anno2+1
			if (index_anno<anno_len-1):
				if (end < contigs[contig1][1]):
					#print("IGR after:" + str(end) + " in " + contig1)
					igr_count += 1
					igr_tuple.append((contig1, (igr_base + str(igr_count)), end + 1, contigs[contig1][1]))
				if (start2 >1):
					#print("IGR before" + str(end) + " in " + contig1)
					igr_count += 1
					igr_tuple.append((contig2, (igr_base + str(igr_count)),1, start2))
			elif index_anno == anno_len-1:
				if (end < contigs[contig1][1]):
					#print("IGR after:" + str(end) + " in " + contig1)
					igr_count += 1
					igr_tuple.append((contig1, (igr_base + str(igr_count)), end + 1, contigs[contig1][1]))
			break
		elif index_anno2==anno_len-1:
			if (start2 <= end + 1):
				end = end2
				index_anno2 += 1
				index_anno=anno_len
				if (end < contigs[contig1][1]):
					#print("IGR after:" + str(end) + " in " + contig1)
					igr_count += 1
					igr_tuple.append((contig1, (igr_base + str(igr_count)), end + 1, contigs[contig1][1]))
			else:
				index_anno = anno_len
				#print("Non-overlap found:" + str(end) + ":" + str(start2))
				igr_count += 1
				igr_tuple.append((contig1, (igr_base + str(igr_count)), end + 1, start2 - 1))
				if (end2 < contigs[contig1][1]):
					#print("IGR after:" + str(end2) + " in " + contig1)
					igr_count += 1
					igr_tuple.append((contig1, (igr_base + str(igr_count)), end2 + 1, contigs[contig1][1]))
					break
for igr in igr_tuple:
	gff3_out.append(igr[0] + "\tmanual\tgene\t" + str(igr[2]) + "\t" + str(igr[3]) + "\t.\t+\t.\tID=" + igr[1] + ";\n")
with open(argstr.output, 'w') as g_out:
	for line in gff3_out:
		g_out.write(line)
