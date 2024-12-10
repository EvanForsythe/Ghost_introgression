#This is a python script for testing ghost vs ingroup introgression

# inport needed modules

import os
import re
import sys
import glob
import shutil
import numpy as np 
import pandas as pd
import subprocess
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
import Bio.Align
from Bio.Align import MultipleSeqAlignment
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator

#Homemade modules
from internal_functions import check_dat

#get the command line arguments
parser = argparse.ArgumentParser(description='Script for performing ghost int testing') 

parser.add_argument('-out', '--outgroup', type=str, metavar='', required=True, help='unique str found in outgroup seqID')
parser.add_argument('-P1', '--pop1', type=str, metavar='', required=True, help='unique str found in pop1 seqID')
parser.add_argument('-P2', '--pop2', type=str, metavar='', required=True, help='unique str found in pop2 seqID')
parser.add_argument('-P3', '--pop3', type=str, metavar='', required=True, help='unique str found in pop3 seqID')
parser.add_argument('-i', '--indir', type=str, metavar='', required=True, help='full path to input directory (should end in "/")')
parser.add_argument('-j', '--job', type=str, metavar='', required=True, help='name for the job to create file names')
parser.add_argument('-t', '--threads', type=int, metavar='', required=True, help='number of threads to use for tree inference')

args = parser.parse_args() 

outgroup = args.outgroup
P1 = args.pop1
P2 = args.pop2
P3 = args.pop3
indir = args.indir
job = args.job
threads = args.threads

# Check the types of variables
check_dat(outgroup, "outgroup", str)
check_dat(P1, "P1", str)
check_dat(P2, "P2", str)
check_dat(P3, "P3", str)
check_dat(indir, "indir", str)
check_dat(job, "job", str)
check_dat(threads, "threads", int)


print(f"outgroup taxon: {outgroup}")
print(f"population 1 taxon: {P1}")
print(f"population 2 taxon: {P2}")
print(f"population 3 taxon: {P3}") 

#create an output directory
out_dir = "OUT_"+job+"/"

#Make a directory for storing stats
if os.path.isdir(out_dir):
    while True:
        user_input = input("This jobname already exists. Would you like to overwrite? (y/n) \n")
        if user_input == 'y':
            print("Clearing contents of " + out_dir + " All output files will be written to this folder\n")
            shutil.rmtree(out_dir)
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

#make a directory for the trees
os.makedirs(out_dir+"alns_and_trees/")

#get a list of seq files to analyze 
all_file_names = glob.glob(indir+"*.fasta")

#create an empty dataframe
node_depth_df = pd.DataFrame()

#Add columns
node_depth_df['gene_name'] = [] 
node_depth_df['node_depth'] = []
node_depth_df['topology'] = []

#counter is for testing purposes 
counter = 0

print("Creating pruned alignments...")

#Loop through each file
for filename in all_file_names:
    #next three lines are for quick testing
    #counter += 1
    #if counter > 100:
    #    break
    sequences = SeqIO.parse(filename, 'fasta')

    #Create an empty list    
    keeper_rec_list=[]

    #loop through each sequence in the file    
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
            
    #Create an alignment of sequences
    keeper_aln = Bio.Align.MultipleSeqAlignment(keeper_rec_list)
    AlignIO.write(keeper_aln, out_dir+"alns_and_trees/pruned_"+filename.replace(indir, ""), "fasta")


#get a list of pruned seq files to analyze 
pruned_file_names = glob.glob(out_dir+"alns_and_trees/"+"pruned*.fasta")

print(f"Inferring gene trees for {len(pruned_file_names)} alignments...")

for file in pruned_file_names:
    tree_command = ["iqtree", "-s", file, "-m", "TEST", "-nt", str(threads)]
    #print(tree_command)
    subprocess.run(tree_command, stdout = subprocess.DEVNULL)

#get a list of tree files to analyze
tree_file_names = glob.glob(out_dir+"alns_and_trees/"+"pruned*.treefile")

for tree_file in tree_file_names:
    temp_tree = Phylo.read(tree_file, "newick")
    temp_tree.root_with_outgroup({"name": "Outgroup"}) 
    all_terminal_branches = temp_tree.get_terminals()
    
    for t in all_terminal_branches:
        if t.name == "P1":
            P1_temp=t
            P1_bl_temp=t.branch_length    
        elif t.name == "P2":
            P2_temp=t
            P2_bl_temp=t.branch_length
        elif t.name == "P3":
            P3_temp=t
            P3_bl_temp=t.branch_length
        else:
            out_temp=t
            
    P1_and_P2=[P1_temp, P2_temp]
    P1_and_P3=[P1_temp, P3_temp]
    P2_and_P3=[P2_temp, P3_temp]
    
    if bool(temp_tree.is_monophyletic(P1_and_P2)):
        topo_str = "12top"
        node_depth=P1_bl_temp + P2_bl_temp
    elif bool(temp_tree.is_monophyletic(P1_and_P3)):
        topo_str = "13top"
        node_depth=P1_bl_temp + P3_bl_temp
    elif bool(temp_tree.is_monophyletic(P2_and_P3)):
        topo_str = "23top"
        node_depth=P2_bl_temp + P3_bl_temp
    else:
        topo_str = "Unknown"
        node_depth="Unknown" 
        
    node_depth_df.loc[len(node_depth_df.index)] = [tree_file, node_depth, topo_str]
    #print(f"topology: {topo_str} node depth: {node_depth}")

print("Writing node depth file...")

node_depth_df.to_csv(out_dir+"node_depths.csv" , index=False)

#node_depth_df = pd.read_csv(out_dir+"node_depths.csv")


# write to a quant file
quant_log_file = "Quant_results_log.tsv"

#Create the quantitative data log file
if not os.path.isfile(quant_log_file):
	with open(quant_log_file, "a") as f:
		f.write("Job_name\tn_12top\tn_23top\tn_13top\tn_unknown_top\td_stat\tavg_12top\tavg_23top\tdelta\trep_pop_mean\tz_score\tp_value\n")


topo_list  = list(node_depth_df["topology"])

#print(topo_list)

#create counter of tree topologies
ticker12 = 0
ticker13 = 0
ticker23 = 0
tickerunknown = 0

for i in topo_list: 
    if i == "12top":
        ticker12 += 1 
    elif i == "13top":
        ticker13 += 1
    elif i == "23top":
        ticker23 += 1
    else:
        tickerunknown += 1 

#Make sure that none of the tickers are
if ticker12 < 5 or ticker23 < 5:
    print("ERROR: one of the required topologies was present less than 5 times")
    with open (quant_log_file, "a") as f:
         f.write(f"{job}\t{ticker12}\t{ticker23}\t{ticker13}\t{tickerunknown}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
    sys.exit()

# Calculate the D-statistic
d_stat = (ticker23-ticker13)/(ticker23+ticker13) 


print(f"Number of 12-top trees: {ticker12}")
print(f"Number of 23-top trees: {ticker23}")
print(f"Number of 13-top trees: {ticker13}")
print(f"D-statistic based on tree counts is: {d_stat}")

###### start generating statistics 

#create violin Plot
#plt.rcdefaults()
sns.violinplot(x = "topology", y = "node_depth", data = node_depth_df.drop(node_depth_df[node_depth_df['topology'] == '13top'].index), split = True)
plt.ylim(0,0.5)
plt.savefig(out_dir+"violin_plot.pdf")
plt.close()

avg_12top=np.average(list(node_depth_df.loc[node_depth_df['topology'] == '12top']["node_depth"]))
avg_23top=np.average(list(node_depth_df.loc[node_depth_df['topology'] == '23top']["node_depth"]))
delta=avg_12top-avg_23top

#Bootstrap resample

print("Running bootstrap resampling for delta statistic")

#create empty list for bootstrap deltas
delta_reps=[]

for i in range(0,100):
    rand_df=node_depth_df.sample(n=len(node_depth_df.index), replace=True)

    rand_avg_12top=np.average(list(rand_df.loc[rand_df['topology'] == '12top']["node_depth"]))
    rand_avg_23top=np.average(list(rand_df.loc[rand_df['topology'] == '23top']["node_depth"]))
    rand_delta=rand_avg_12top-rand_avg_23top
    delta_reps.append(rand_delta)

#conduct z_test and get a p_value

# Given information
sample_mean = 0
population_mean = np.average(delta_reps)
print('pop mean:',population_mean)

#population_std = stats.tstd(delta_reps)
population_std = np.std(delta_reps, ddof=0)

print('pop std:',population_std)
sample_size = len(delta_reps)
print('sample size:',sample_size)

# compute the z-score
z_score = (sample_mean-population_mean)/(population_std/np.sqrt(sample_size))
print('Z-Score :',z_score)


# P-Value : Probability of getting less than a Z-score
p_value = (1-stats.norm.cdf(abs(z_score)))*2   #two-tailed test
print(f'p-value :{p_value}')

#create kernal density plot
plt.rcdefaults()
sns.kdeplot(delta_reps)
plt.axvline(x = 0, color = 'b')
plt.title(f"Distribution mean: {population_mean}\np_value: {p_value}")
plt.savefig(out_dir+"delta_plot.pdf")
plt.close()


#Write to the quant file
with open (quant_log_file, "a") as f:
	f.write(f"{job}\t{ticker12}\t{ticker23}\t{ticker13}\t{tickerunknown}\t{d_stat}\t{avg_12top}\t{avg_23top}\t{delta}\t{population_mean}\t{z_score}\t{p_value}\n")


print("Done!")

