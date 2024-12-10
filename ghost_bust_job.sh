#!/bin/bash
#SBATCH --job-name=prop05
#SBATCH --ntasks-per-node=32
#SBATCH --time=200:0:0
#SBATCH --output=prop05.out
#SBATCH --error=prop05.err
#SBATCH --mail-user=evan.forsythe@osucascades.edu
#SBATCH --mail-type=END

#conda activate ghost_int
#sbatch -p forsythe.q -A forsythe ghost_bust_job.sh

#python Data_simulations.py -j prop01 -p 0.4 --ghost
#python Data_simulations.py -j prop02 -p 0.2 --ghost
#python Data_simulations.py -j prop03 -p 0.3 --ghost
#python Data_simulations.py -j prop04 -p 0.4 --ghost
#python Data_simulations.py -j prop05 -p 0.5 --ghost

#python Ghost_buster.py -i OUT_prop01/single_gene_alns/ -j bust_prop01 -t 24 -P1 Pop1 -P2 Pop2 -P3 Pop3 -out Outgroup
#python Ghost_buster.py -i OUT_prop02/single_gene_alns/ -j bust_prop02 -t 24 -P1 Pop1 -P2 Pop2 -P3 Pop3 -out Outgroup
#python Ghost_buster.py -i OUT_prop03/single_gene_alns/ -j bust_prop03 -t 24 -P1 Pop1 -P2 Pop2 -P3 Pop3 -out Outgroup
#python Ghost_buster.py -i OUT_prop04/single_gene_alns/ -j bust_prop04 -t 24 -P1 Pop1 -P2 Pop2 -P3 Pop3 -out Outgroup
python Ghost_buster.py -i OUT_prop05/single_gene_alns/ -j bust_prop05 -t 24 -P1 Pop1 -P2 Pop2 -P3 Pop3 -out Outgroup
