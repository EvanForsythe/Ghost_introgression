#This is a python script for testing ghost vs ingroup introgression



# inport needed modules

import os
import glob
#import numpy as np 
#import pandas as pd
from Bio import SeqIO
from Bio import AlignIO



###get a list of seq files to analyze 
	#store as python list object

file_name_str = "Examp_files/*.phy"
#print(type(file_name_str))

all_file_names = glob.glob(file_name_str)

###convert each file to FASTA format and write a new version of the file. **WORKING**

#### NOTE: Everything below this is under constructions
'''
#print(type(all_file_names))

###print the number of files in this list.

#print(len(all_file_names))

###read in the first file in the list and store as a python alignment object.

name_0=all_file_names[0]
print("our first file is named "+str(name_0)+" and is stored as a "+str(type(name_0))+" object.")

align = AlignIO.read(name_0 , "phylip")
print(type(align))
print(align)

###print the number of sites in the alignment. 

'''
