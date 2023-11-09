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


'''
#Compare the two lists so there's not redundancy.
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
import Bio.Align
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
'''

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
'''
#THIS BLOCK MAY NEED TO BE DELETED

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
'''


#### Below here is the code we developed in a jupyter notebook.


#Make adjustments to the following lines
'''
filename = "/home/forsythe/AlignIO_testing/At_1G06920.fasta" 
sequences = SeqIO.parse(filename, 'fasta')

P1 = "Bs_"
P2 = "Cr_"
P3 = "At_"
outgroup = "Es_"
'''

keeper_rec_list=[]

for record in sequences:
    if P1 in record.id:
        #print(record.id)
        P1_rec_temp = record
        P1_rec_temp.id = "P1"
        keeper_rec_list.append(P1_rec_temp)
    elif P2 in record.id:
        #print(record.id)
        P2_rec_temp = record
        P2_rec_temp.id = "P2"
        keeper_rec_list.append(P2_rec_temp)
    elif P3 in record.id:
        #print(record.id)
        P3_rec_temp = record
        P3_rec_temp.id = "P3"
        keeper_rec_list.append(P3_rec_temp)
    elif outgroup in record.id:
        #print(record.id)
        out_rec_temp = record
        out_rec_temp.id = "Outgroup"
        keeper_rec_list.append(out_rec_temp)

keeper_aln = Bio.Align.MultipleSeqAlignment(keeper_rec_list)

calculator = DistanceCalculator('identity')

distance_matrix = calculator.get_distance(keeper_aln)
print(distance_matrix)

#Create a weird constructor object
constructor = DistanceTreeConstructor()

# Construct the phlyogenetic tree using NJ algorithm
NJTree = constructor.nj(distance_matrix)
# Draw the phlyogenetic tree using terminal
Bio.Phylo.draw_ascii(NJTree)

all_terminal_branches = NJTree.get_terminals()

P1_and_P2=[]
P1_and_P3=[]
P2_and_P3=[]

for t in all_terminal_branches:
    if t.name == "P1":
        P1_temp=t
    elif t.name == "P2":
        P2_temp=t
    elif t.name == "P3":
        P3_temp=t
    
    P1_and_P2=[P1_temp, P2_temp]
    P1_and_P3=[P1_temp, P3_temp]
    P2_and_P3=[P2_temp, P3_temp]

if bool(NJTree.is_monophyletic(P1_and_P2)):
    topo_str = "12top"
elif bool(NJTree.is_monophyletic(P1_and_P3)):
    topo_str = "13top"
elif bool(NJTree.is_monophyletic(P2_and_P3)):
    topo_str = "23top"
else:
    topo_str = "Unknown"

print(topo_str)
