#!/Users/shounak/anaconda3/bin/python3.7

# This program prepares an input for OPaver, given a CSV file with differential methylation data
# Opaver input file is a TSV with four columns
#Annotation string, KO numbers separated by ';',meth_diff,p-value
import argparse
#function definition start
def collect_KO(genome_ID,ko_file,cut_off): #This function collects KO numbers for a given genome ID/accession number; shortlist created based on E-value cutoff
	result=""
	genome_ID_trimmed=genome_ID.replace("* ","")
	for an_ID in genome_ID_trimmed.split(' '):
		an_ID_trim=an_ID.split('-')[0]
		try:
			ko_f=open(ko_file,'r') #try to open the file specified by ko_file (an absolute path)
		except FileNotFoundError:
			print ("Error opening KO file:",ko_file) #if something goes wrong; inform the user and
			quit() #quit the program
	#now that the KO file has been opened, we will collect KO numbers from the master list if their E-values fall below specified cutoff
		ko_list=ko_f.readlines()
		for a_line in ko_list:
			field_a=a_line.split(',') #split the comma separated fields
			KO_list_genome_ID=field_a[0]#.split('.')[0]
			if (an_ID_trim==KO_list_genome_ID):
				if (field_a[4]!='-'):
					if (float(field_a[4])<=cut_off):
						result=result+field_a[1]+";" # collect the KO numbers, separated by a ';' for OPaver
		ko_f.close()
	return result
#function definition end

#Begin main program
#We will set up the following command line arguments
#-f --csvfile : path to CSV input file 
#-n --id_column : column to look in for accession IDs
#-k --kofile: path to CSV file containing KO numbers for each genome_ID
#-e --eval_lim: Kofam HMM E-value cutoff
#-o --output : path to CSV output file

#Let's see what the user has given us

optParse=argparse.ArgumentParser() #Create an argument parsing "object" called optParse
optParse.add_argument("-f","--csvfile",help="path to CSV input file")
optParse.add_argument("-n","--id_column",help="column to look in for accession IDs")
optParse.add_argument("-k","--kofile",help="path to CSV file containing KO numbers for each genome_ID")
optParse.add_argument("-e","--eval_lim",help="Kofam HMM E-value cutoff")
optParse.add_argument("-o","--output",help="path to CSV output file")

argstr=optParse.parse_args() # store the arguments given/passed by the user
def_cutoff=1e-12
if argstr.eval_lim is not None:
	try:
		def_cutoff=float(argstr.eval_lim)
	except AttributeError:
		print ("Please check the cut off E value supplied")
		quit()
#Now, get to work and try to open the files
try:
	csv_f=open(argstr.csvfile,'r') # open the input CSV file in read-only mode
except FileNotFoundError:
	print ("Error opening input CSV file:",argstr.csvfile) #if something goes wrong; inform the user and
	quit() #quit the program
	
out_f=open(argstr.output,'w')# open the output CSV file in write-only mode

#Now let's start reading the csv file, line by line
# We have three types of inputs: promoters, exons and introns as ends of filename
# Based on these types, we need to look in the right column for accession IDs (in this case, NSC_XXXXXX)
# Some rows will have multiple accession IDs
id_col=int(argstr.id_column)-1 # assign variable to store the column number to look up (and convert it to python 'list numbering''; basically human column value minus 1
csv_line_ct=1
for line in csv_f:
	#store field values in an array
	inp_field=line.split(',')
	print(len(inp_field))
	if (inp_field[id_col] and inp_field[id_col]!="\n" and inp_field[id_col]!=" "):
		inp_field[id_col]=inp_field[id_col].replace('"','')
		if csv_line_ct!=1:
			sub_inp=inp_field[id_col].rstrip("\n").rstrip("_P").rstrip("_T")# keeping this from old script so we can handle transcript level input later
			#annotation=inp_field[id_col]+"-()-"+inp_field[0]+"-()-"+inp_field[1]+"-()-"+inp_field[2]+"-()-"
			annotation = inp_field[59]
			ko_str=collect_KO(sub_inp,argstr.kofile,def_cutoff)
			#print(annotation+"\t"+ko_str+"\t"+inp_field[21]+"\t"+inp_field[22]+"\n")
			#out_f.write(annotation+"\t"+ko_str+"\t"+inp_field[21]+"\t"+inp_field[22]+"\n")
			print(annotation + "\t" + ko_str + "\t" + inp_field[56] + "\t" + "0\n")
			out_f.write(annotation + "\t" + ko_str + "\t" + inp_field[56] + "\t" + "0\n")
		else:
			out_f.write("Annotation\tKO\tmeth_diff\tpvalue\n")
	csv_line_ct=csv_line_ct+1
csv_f.close()
out_f.close()
