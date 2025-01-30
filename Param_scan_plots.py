#Load packages
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#Import data into a dataframe
sim_df = pd.read_table('/scratch/pappab/Ghost_project/Ghost_introgression/Sim_stats_log.tsv')
quant_df = pd.read_table('/scratch/pappab/Ghost_project/Ghost_introgression/Quant_results_log.tsv')


scaling_quant_df = quant_df[quant_df["Job_name"].str.contains("scale_")]


sf_list=list()
rep_list=list()
Int_type_list=list()

for row in scaling_quant_df.index:
    temp_jobname = scaling_quant_df.loc[row]["Job_name"]
    temp_sf=temp_jobname.split("_")[4]
    temp_rep=temp_jobname.split("_")[-1]
    temp_type=temp_jobname.split("_")[1]
    
    sf_list.append(float(temp_sf))
    rep_list.append(int(temp_rep))
    Int_type_list.append(str(temp_type))
    #scaling_quant_df.loc[row]["Scaling_factor"] = temp_sf
    #scaling_quant_df.loc[row]["Replicate"] = temp_rep

scaling_quant_df["Scaling_factor"] = sf_list
scaling_quant_df["Replicate"] = rep_list
scaling_quant_df["Int_type"] = Int_type_list



#Create plot of scaling factor
sns.scatterplot(
    data=scaling_quant_df, 
    x="Scaling_factor", 
    y="rep_pop_mean", 
    hue="Int_type",  # Categorical column for color mapping
    palette="viridis",  # Use the same colormap as before
    edgecolor="black"  # Optional for better visibility
)

# Add legend, labels, and title
plt.legend(title='Introgression type')
plt.xlabel('Scaling Factor')
plt.ylabel('Delta (average from bootstrap resampling)')
plt.title('Scaling factor parameter scan plot')

# Save the figure
plt.savefig('Scaling_factor_fig.pdf', format="pdf")

plt.show()
plt.close()

#Introgression time plot

#create Data frame
int_time_quant_df = quant_df[quant_df["Job_name"].str.contains("int_time_")]

#Add columns to data frame
int_time_list=list()
rep_list=list()
Int_type_list=list()

for row in int_time_quant_df.index:
    temp_jobname = int_time_quant_df.loc[row]["Job_name"]
    temp_time=temp_jobname.split("_")[5]
    temp_rep=temp_jobname.split("_")[-1]
    temp_type=temp_jobname.split("_")[1]
    
    int_time_list.append(float(temp_time))
    rep_list.append(int(temp_rep))
    Int_type_list.append(str(temp_type))
    #scaling_quant_df.loc[row]["Scaling_factor"] = temp_sf
    #scaling_quant_df.loc[row]["Replicate"] = temp_rep

int_time_quant_df["Introgresion_time"] = int_time_list
int_time_quant_df["Replicate"] = rep_list
int_time_quant_df["Int_type"] = Int_type_list


#Create plot of results
sns.scatterplot(
    data=int_time_quant_df, 
    x=("Introgresion_time"), 
    y="rep_pop_mean", 
    hue="Int_type",  # Categorical column for color mapping
    palette="viridis",  # Use the same colormap as before
    edgecolor="black",  # Optional for better visibility
    #x_jitter=1
    )


# Add legend, labels, and title
#plt.xlim(4000, 80000)
plt.legend(title='Introgression type')
plt.xlabel("Introgresion time")
plt.ylabel('Delta (average from bootstrap resampling)')
plt.title('Introgression time parameter scan plot')

# Save the figure
plt.savefig('Introgression_time_fig.pdf', format="pdf")

plt.show()
plt.close()
