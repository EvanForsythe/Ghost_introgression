#!/bin/bash
#SBATCH --job-name=python_array           # Job name
#SBATCH --output=output_%A_%a.out         # Output file name
#SBATCH --error=error_%A_%a.err           # Error file name
#SBATCH --time=24:00:00                   # Time limit (hh:mm:ss)
#SBATCH --ntasks=1                        # Number of tasks (usually 1 for array jobs)
#SBATCH --cpus-per-task=2                 # Number of CPU cores per task
#SBATCH --mem=4G                          # Memory per node
#SBATCH --array=1-40                      # Array range 


#conda activate ghost_int
#sbatch -p forsythe.q -A forsythe ghost_array_2D.sh


# Define the number of unique values for r and m
r_values_count=5
s_values_count=8

# Calculate r_index and m_index based on SLURM_ARRAY_TASK_ID
task_id=$SLURM_ARRAY_TASK_ID
s_index=$(( (task_id - 1) / r_values_count + 1 ))
r_index=$(( (task_id - 1) % r_values_count + 1 ))

# Generate scaling factor
scaling_factor=$(awk "BEGIN {print $s_index * 0.25}")

#echo "scaling_factor: $scaling_factor"

#Set the default divergences
default_tI=40000
default_t2=80000
default_t3=120000
default_tG=160000
default_t4=200000

#Set the scaled divergences
# Scale each divergence
scaled_tI=$(awk "BEGIN {print $default_tI * $scaling_factor}")
scaled_t2=$(awk "BEGIN {print $default_t2 * $scaling_factor}")
scaled_t3=$(awk "BEGIN {print $default_t3 * $scaling_factor}")
scaled_tG=$(awk "BEGIN {print $default_tG * $scaling_factor}")
scaled_t4=$(awk "BEGIN {print $default_t4 * $scaling_factor}")

# Output the scaled values
echo "Scaling Factor: $scaling_factor"
echo "Scaled tI: $scaled_tI"
echo "Scaled t2: $scaled_t2"
echo "Scaled t3: $scaled_t3"
echo "Scaled tG: $scaled_tG"
echo "Scaled t4: $scaled_t4"

#Create the jobname ## NOTE CHANGE THE TRUE/GHOST HERE
j_value="true__sf_${scaling_factor}__rep_${r_index}"  # Create j_value string
echo "j_value: $j_value"


# Run the simulation job
echo Executing: python Data_simulations.py -j $j_value -i $scaled_tI -2 $scaled_t2 -3 $scaled_t3 -G $scaled_tG -4 $scaled_t4
# Uncomment below to run the simulation
python Data_simulations.py -j "$j_value" -i "$scaled_tI" -2 "$scaled_t2" -3 "$scaled_t3" -G "$scaled_tG" -4 "$scaled_t4"

# Run the ghost buster job
echo Executing: python Ghost_buster.py -j bust_${j_value} -i OUT_${j_value}/single_gene_alns/ -t 2 -P1 Pop1 -P2 Pop2 -P3 Pop3 -out Outgroup
# Uncomment below to run the ghost buster
python Ghost_buster.py -j "bust_${j_value}" -i "OUT_${j_value}/single_gene_alns/" -t 2 -P1 Pop1 -P2 Pop2 -P3 Pop3 -out Outgroup