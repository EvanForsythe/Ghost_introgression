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
import itertools

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
parser.add_argument('-r', '--Recomb_rate', type=float, metavar='', required=False,default=0.000000001, help='Specify the recomb rate (default = 0.000000001)')
parser.add_argument('-n', '--Ne', type=int, metavar='', required=False, default=10000, help='Specify the effective pop size (Ne) (default = 10000)')
parser.add_argument('-g','--ghost', action='store_true', required=False, help='Add this flag to to simulate ghost introgression. Otherwise introgression will be "true"/"ingroup" introgression from P3 to P2"')


#Time arguments
parser.add_argument('-i', '--t_int', type=int, metavar='', required=False, default=40000 , help='Time of introgression (years ago) (default = 40000)')
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
check_dat(JOBname, "JOBname", str)
check_dat(Seq_len, "Seq_len", int)
check_dat(Prop_int, "Prop_int", float)
check_dat(Mut_rate, "Mut_rate", float)
check_dat(Recomb_rate, "Recomb_rate", float)
check_dat(Ne, "Ne", int)
check_dat(ghost, "ghost", bool)
check_dat(t_int, "t_int", int)
check_dat(t_sp12, "t_sp12", int)
check_dat(t_sp123, "t_sp123", int)
check_dat(t_sp123G, "t_sp123G", int)
check_dat(t_sp123G4, "t_sp123G4", int)

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

# Add introgression and speciation events

# Track introgression events manually
introgression_events = []

#Ask if we're in ghost mode or not
if not ghost:
    demography.add_mass_migration(
        time=t_int, source="Pop2", dest="Pop3", proportion=Prop_int)
    introgression_events.append({"time": t_int, "source": "Pop2", "dest": "Pop3", "proportion": Prop_int})
else:
    demography.add_mass_migration(
        time=t_int, source="Pop1", dest="Ghost", proportion=Prop_int)
    introgression_events.append({"time": t_int, "source": "Pop1", "dest": "Ghost", "proportion": Prop_int})

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

# Count the number of tracts (haplotype blocks) and calculate their average length
total_length = 0
num_tracts = 0

for migration in ts.migrations():
    tract_length = migration.right - migration.left  # Calculate tract length
    total_length += tract_length
    num_tracts += 1

# Calculate the average length
average_length = total_length / num_tracts if num_tracts > 0 else 0

print(f"Number of tracts (haplotype blocks) simulated: {num_tracts}")
print(f"Average tract length (in nucleotides): {average_length:.2f}")
print(f"total genome size is {Seq_len}, meaning each window will be {Seq_len/1000}")

# Generate mutations on the tree sequence
ts_mutes = msprime.sim_mutations(ts, rate=Mut_rate, random_seed=None)

# Calculate total sequence divergence (average pairwise divergence)
def calculate_total_divergence(tree_sequence):
    divergence = 0
    pair_count = 0

    # Iterate over pairs of samples
    for i, j in itertools.combinations(range(tree_sequence.num_samples), 2):
        divergence += tree_sequence.divergence([[i], [j]])  # Pairwise divergence
        pair_count += 1

    average_divergence = divergence / pair_count if pair_count > 0 else 0
    return average_divergence


total_divergence = calculate_total_divergence(ts_mutes)
print(f"Total sequence divergence (average pairwise divergence): {total_divergence:.6f} substitutions/site")

#Get the 'diagnostic SNPs' that indicate introgression
def count_introgressed_mutations(ts, introgression_events):
    """
    Counts the number of mutations that occurred along branches in introgressed regions,
    calculates the average length of introgressed tracts, and computes the total
    length of non-overlapping introgressed tracts.

    Args:
        ts (TreeSequence): The simulated tree sequence.
        introgression_events (list): List of tracked introgression events.

    Returns:
        tuple: A tuple containing:
            - int: Total number of introgressed tracts.
            - float: Average length of introgressed tracts.
            - int: Total number of "diagnostic" mutations in introgressed tracts.
            - float: Total length of non-overlapping introgressed tracts.
    """
    # Generate introgressed tracts based on introgression events
    introgressed_tracts = []
    for migration in ts.migrations():
        for event in introgression_events:
            source_id = next((pop.id for pop in ts.populations() if pop.metadata.get("name") == event["source"]), None)
            dest_id = next((pop.id for pop in ts.populations() if pop.metadata.get("name") == event["dest"]), None)
            if source_id is not None and dest_id is not None:
                if (event["time"] - 1 <= migration.time <= event["time"] + 1 and
                        migration.source == source_id and migration.dest == dest_id):
                    introgressed_tracts.append((migration.left, migration.right))

    # Sort and merge overlapping tracts
    introgressed_tracts.sort()
    merged_tracts = []
    if introgressed_tracts:
        current_start, current_end = introgressed_tracts[0]
        for start, end in introgressed_tracts[1:]:
            if start <= current_end:  # Overlap
                current_end = max(current_end, end)
            else:  # No overlap
                merged_tracts.append((current_start, current_end))
                current_start, current_end = start, end
        merged_tracts.append((current_start, current_end))  # Add the last tract

    # Calculate the total length of introgressed tracts
    total_tract_length = sum(end - start for start, end in merged_tracts)

    # Calculate the average tract length
    num_tracts = len(merged_tracts)
    average_tract_length = total_tract_length / num_tracts if num_tracts > 0 else 0

    # Loop through mutations and check if they are in introgressed regions
    diagnostic_mutations = 0
    for mutation in ts.mutations():
        site_position = ts.site(mutation.site).position
        for left, right in merged_tracts:
            if left <= site_position < right:
                diagnostic_mutations += 1
                break

    return num_tracts, average_tract_length, diagnostic_mutations, total_tract_length


#Call the function
num_int_tracts, average_int_tract_length, num_diagnostic_mutations, total_introgressed_length = count_introgressed_mutations(ts_mutes, introgression_events)

print(f"Number of introgressed tracts: {num_int_tracts}")
print(f"Average length of int tract: {average_int_tract_length}")
print(f"Number of 'diagnostic' mutations in introgressed tracts: {num_diagnostic_mutations}")
print(f"Total length of non-overlapping introgressed tracts: {total_introgressed_length}")

# write to a quant file
quant_log_file = "Sim_stats_log.tsv"

#Create the quantitative data log file
if not os.path.isfile(quant_log_file):
	with open(quant_log_file, "a") as f:
		f.write("JOBname\tSeq_len\tProp_int\tMut_rate\tRecomb_rate\tnum_tracts\taverage_length\ttotal_divergence\tnum_int_tracts\taverage_int_tract_length\ttotal_introgressed_length\tnum_diagnostic_mutations\n")

#Write to the quant file
with open (quant_log_file, "a") as f:
	f.write(f"{JOBname}\t{Seq_len}\t{Prop_int}\t{Mut_rate}\t{Recomb_rate}\t{num_tracts}\t{average_length}\t{total_divergence}\t{num_int_tracts}\t{average_int_tract_length}\t{total_introgressed_length}\t{num_diagnostic_mutations}\n")

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
    print("Done!")

if __name__ == "__main__":
    sliding_window_fasta(out_dir+"full_chromosome_renamed.fa", out_dir+"single_gene_alns/")
