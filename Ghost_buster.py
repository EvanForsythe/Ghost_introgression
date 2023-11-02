#This is a python script for testing ghost vs ingroup introgression



# inport needed modules

import os
import re
import glob
#import numpy as np 
#import pandas as pd
import argparse
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

#get the command line arguments

parser = argparse.ArgumentParser(description='Script for performing ghost int testing') 

parser.add_argument('-out', '--outgroup', type=str, metavar='', required=True, help='unique str found in outgroup seqID')
parser.add_argument('-P1', '--pop1', type=str, metavar='', required=True, help='unique str found in pop1 seqID')
parser.add_argument('-P2', '--pop2', type=str, metavar='', required=True, help='unique str found in pop2 seqID')
parser.add_argument('-P3', '--pop3', type=str, metavar='', required=True, help='unique str found in pop3 seqID')

args = parser.parse_args() 

outgroup = args.outgroup
pop1 = args.pop1
pop2 = args.pop2
pop3 = args.pop3

print(outgroup)
print(pop1)
print(pop2)
print(pop3)

#get a list of seq files to analyze 

all_file_names = glob.glob("Examp_files/*.fasta")
#print(all_file_names)

###read in the first file in the list and store as a python alignment object.

my_file = all_file_names[0]
#print(my_file)

###create fasta file from phy file

alignment = AlignIO.read(my_file, "fasta")
#print(type(alignment)
print(type(alignment[0]))

###create an empty list

record_ids = []

#pruned_alignment = MultipleSeqAlignment()

###add seqrecord id to the empty list

for i in alignment:
	record_ids.append(i.id)
	if outgroup in i.id:
		outgroup_id = i.id
		#pruned_alignment = MultipleSeqAlignment(i)
	elif pop1 in i.id:
		pop1_id = i.id
		#pruned_alignment.append(i)
	elif pop2 in i.id:
		pop2_id = i.id
		#pruned_alignment.append(i)
	elif pop3 in i.id:
		pop3_id = i.id
		#pruned_alignment.append(i)
	else: 
		print(f"trimming {i.id} from the alignment")

#print(pruned_alignment)




#Qcount = SeqIO.write(records, my_file.replace('.phy',''), "fasta")

###print the number of sites in the alignment. 


