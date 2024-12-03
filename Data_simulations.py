#!/usr/bin/env python3

# script for running introgression simulations
import sys
import random
import collections
import subprocess
import matplotlib.pyplot as plt
import msprime
import numpy as np
import dataclasses
import tskit
import os
import re
import argparse
import math
import pandas as pd
import shutil as sh
from Bio import SeqIO

#Homemade modules
from internal_functions import check_dat

#Set wd
working_dir = sys.path[0]+'/' 
os.chdir(working_dir)

#Set up an argumanet parser
parser = argparse.ArgumentParser(description='Script for simulating introgression')

#Add arguments
parser.add_argument('-j', '--JOBname', type=str, metavar='', required=True, help='Unique job name for this run of this script. Avoid including spaces or special characters ("_" is ok)') 
parser.add_argument('-s', '--Seq_len', type=int, metavar='', required=False, default=10000000, help='Specify an interger to set length of total simulateed alignment (default = 10000000')
parser.add_argument('-p', '--Prop_int', type=float, metavar='', required=False, default=0.2, help='Specify the proportion of genome to be introgressed with each introgression event (default = 0.2)')
parser.add_argument('-m', '--Mut_rate', type=float, metavar='', required=False, default=0.0000001, help='Specify the mutation rate (default = 0.0000001)')
parser.add_argument('-r', '--Recomb_rate', type=float, metavar='', required=False,default=0.0000000001, help='Specify the recomb rate (default = 0.0000000001)')
parser.add_argument('-n', '--Ne', type=int, metavar='', required=False, default=10000, help='Specify the effective pop size (Ne) (default = 10000)')
parser.add_argument('-g','--ghost', action='store_true', required=False, help='Add this flag to to simulate ghost introgression. Otherwise introgression will be "true"/"ingroup" introgression from P3 to P2"')


#Time arguments
parser.add_argument('-1', '--t_int', type=int, metavar='', required=False, default=40000 , help='Time of introgression (years ago) (default = 40000)')
parser.add_argument('-2', '--t_sp12', type=int, metavar='', required=False,default=80000 , help='Time of first most recent speciation (years ago) (default = 80000)')
parser.add_argument('-3', '--t_sp123', type=int, metavar='', required=False,default=120000 , help='Time of second most recent speciation (default = 120000)')
parser.add_argument('-G', '--t_sp123G', type=int, metavar='', required=False, default=160000 , help='Time of third  most recent speciation (default = 160000)')
parser.add_argument('-4', '--t_sp123G4', type=int, metavar='', required=False, default=200000 , help='Time of forth  most recent speciation (default = 200000)')

#Define the parser
args = parser.parse_args()

JOBname=args.JOBname
Seq_len=args.Seq_len
Prop_int=args.Prop_int
Mut_rate=args.Mut_rate
Recomb_rate=args.Recomb_rate
Ne=args.Ne
ghost=args.ghost
t_int=args.t_int
t_sp12=args.t_sp12
t_sp123=args.t_sp123
t_sp123G=args.t_sp123G
t_sp123G4=args.t_sp123G4


# Apply the function to check status of all user inputs
check_dat("JOBname", str)
check_dat("Seq_len", int)
check_dat("Prop_int", float)
check_dat("Mut_rate", float)
check_dat("Recomb_rate", float)
check_dat("Ne", int)
check_dat("ghost", bool)
check_dat("t_int", int)
check_dat("t_sp12", int)
check_dat("t_sp123", int)
check_dat("t_sp123G", int)
check_dat("t_sp123G4", int)

#Store output dir as a variable
out_dir= 'OUT_'+JOBname+'/'

#Create the output folder
#Make a directory for storing stats
if os.path.isdir(out_dir):
    while True:
        user_input = input("This jobname already exists. Would you like to overwrite? (y/n) \n")
        if user_input == 'y':
            print("Clearing contents of " + out_dir + " All output files will be written to this folder\n") 
            sh.rmtree(out_dir)
            os.makedirs(out_dir)
            break
        if user_input == 'n':
            print('Unique jobname required. Exiting...')
            sys.exit()
            break
        else:
            print("Command not recognized.\n")     
else:
    os.makedirs(out_dir) 
    print('created folder: '+out_dir+'\nAll output files will be written to this folder\n')


#Get list of taxa
taxa_names=["Pop1", "Pop2", "Pop3", "Ghost", "Outgroup"]

#set up a demographic history 

#Setup the simulations
#time_units = 1000 / 25  # Conversion factor for kya to generations
demography = msprime.Demography()


#Loop through taxa and add each as a population
for t_name in taxa_names:
    demography.add_population(name=t_name, initial_size=Ne)
    print(f'Adding population: {t_name}')


#t_int=40000
#t_sp12=80000
#t_sp123=120000
#t_sp1234=200000

#Ask if we're in ghost mode or not
if not ghost:
    demography.add_mass_migration(
        time=t_int, source="Pop2", dest="Pop3", proportion=Prop_int)
else:
    demography.add_mass_migration(
        time=t_int, source="Pop2", dest="Ghost", proportion=Prop_int)

# Speciation event
demography.add_mass_migration(
    time=t_sp12, source="Pop2", dest="Pop1", proportion=1)

# Speciation event
demography.add_mass_migration(
    time=t_sp123, source="Pop3", dest="Pop1", proportion=1)

# Speciation event
demography.add_mass_migration(
    time=t_sp123G, source="Ghost", dest="Pop1", proportion=1)

# Speciation event
demography.add_mass_migration(
    time=t_sp123G4, source="Pop1", dest="Outgroup", proportion=1)



ts = msprime.sim_ancestry(
    recombination_rate=Recomb_rate,
    sequence_length=Seq_len, 
    samples=[
        msprime.SampleSet(1, ploidy=1, population="Pop1"),
        msprime.SampleSet(1, ploidy=1, population="Pop2"),
        msprime.SampleSet(1, ploidy=1, population="Pop3"),
        msprime.SampleSet(1, ploidy=1, population="Ghost"),
        msprime.SampleSet(1, ploidy=1, population="Outgroup"),
    ],
    demography = demography,
    record_migrations=True,  # Needed for tracking segments.
)


# Generate mutations on the tree sequence
ts_mutes = msprime.sim_mutations(ts, rate=Mut_rate, random_seed=None)


#write a fasta file
ts_mutes.write_fasta(out_dir+"full_chromosome.fa", reference_sequence=tskit.random_nucleotides(ts.sequence_length))

### Pseudo code for editing file:
#Create file handle for the fasta file that we wrote (open for reading)

#Create file handle for a new file (open for 'appending')



# Define the mapping from old IDs to new IDs
id_mapping = {
    'n0': 'Pop1',
    'n1': 'Pop2',
    'n2': 'Pop3',
    'n3': 'Ghost',
    'n4': 'Outgroup'
}

# Read the input FASTA file and write to a new file with updated IDs
with open(out_dir+"full_chromosome.fa", 'r') as infile, open(out_dir+"full_chromosome_renamed.fa", 'w') as outfile:
    for record in SeqIO.parse(infile, 'fasta'):
        # Replace the sequence ID based on the mapping
        if record.id in id_mapping:
            record.id = id_mapping[record.id]
            record.description = ''  # Optionally clear the description
        # Write the updated record to the output file
        SeqIO.write(record, outfile, 'fasta')

# Use a sliding window to create single-gene alignments
os.makedirs(out_dir+"single_gene_alns/")

def sliding_window_fasta(input_fasta, output_dir):
    """
    Creates 1000 smaller alignment files using a non-overlapping sliding window approach.

    Args:
    input_fasta (str): Path to the input FASTA file containing long sequence alignments.
    output_dir (str): Directory to store the output FASTA files.

    Returns:
    None
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Read the input FASTA file and get all sequences
    records = list(SeqIO.parse(input_fasta, "fasta"))
    
    if len(records) == 0:
        print("No sequences found in the input FASTA file.")
        return
    
    # Determine the length of the alignment (assuming all sequences have the same length)
    alignment_length = len(records[0].seq)
    
    # Calculate window length as 1/1000 of the alignment length
    win_len = math.ceil(alignment_length / 1000)

    # Iterate over the alignment in windows
    for window_num in range(1000):
        start = window_num * win_len
        end = min(start + win_len, alignment_length)  # Ensure the last window doesn't exceed alignment length
        
        # Create a new list to store windowed sequences
        window_records = []

        # Extract the windowed sequences for each record
        for record in records:
            window_seq = record.seq[start:end]
            window_record = record[start:end]
            window_records.append(window_record)

        # Write the windowed sequences to a new FASTA file
        output_file = os.path.join(output_dir, f"window_{window_num + 1:04}.fasta")
        SeqIO.write(window_records, output_file, "fasta")

        # Stop if end has reached the alignment length
        if end >= alignment_length:
            break


if __name__ == "__main__":
    sliding_window_fasta(out_dir+"full_chromosome_renamed.fa", out_dir+"single_gene_alns/")