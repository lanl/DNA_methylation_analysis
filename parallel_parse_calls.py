#!/bin/python3
# This program figures out the number of CPU's available and farms out
# Opaver tasks to persistent process handles to each CPU
# Fortunately there is a multiprocessing module in Python that allows us to
# queue tasks to be handed to each CPU
# We will set the maximum size of the queue(JoinableQueue) to use_cores
#
import argparse 
import multiprocessing
import os
from multiprocessing import Process,JoinableQueue,Pool
import subprocess #To actually run command line programs from this Python script
#Let us define a subroutine that runs the OPaver tasks
#and call it from the main program that runs multiple instances of the
#subroutine
def run_parsing_tasks(file_str): # Function to run OPaver.pl and OPaver.r
#start
	print("Working on "+file_str)
	os.system(file_str)
#end
if (__name__=='__main__'):
#Read list of TSV files to process from metadata file
	optParse=argparse.ArgumentParser() #Create an argument parsing "object" called optParse
	optParse.add_argument("-f","--metafile",help="path to metadata file")
	argstr=optParse.parse_args()
#Let's process the metadata file now
	meta_file=open(argstr.metafile,'r')
#Collect all the lines
	meta_lines=meta_file.readlines()
#Close the metadata file
	meta_file.close()
# and process them
	#with JoinableQueue(int(multiprocessing.cpu_count()/2)) as Q1:
	#q1=JoinableQueue(int(multiprocessing.cpu_count()/2))
	p1=Pool(int(multiprocessing.cpu_count()/4))
	p1.map(run_parsing_tasks,meta_lines)
	#for task in list_items:
		#q1.put(run_Opaver_tasks(task),True)
	print ("We used "+str(int(multiprocessing.cpu_count()/4))+" CPUs for this")
#We are all done
	print("Finished processing all files!")

