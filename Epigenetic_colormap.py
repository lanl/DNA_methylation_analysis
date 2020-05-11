import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import argparse
# sphinx_gallery_thumbnail_number = 2
# Read the heatmap data into pandas dataframe
optParse=argparse.ArgumentParser() #Create an argument parsing "object" called optParse
optParse.add_argument("-f","--datafile",help="path to CSV input file")
optParse.add_argument("-m","--map",help="color map to use, see matplotlib.cmap for info")
optParse.add_argument("-o","--output",help="path to heatmap output")
argstr=optParse.parse_args()
df=pd.read_csv(argstr.datafile,sep=',',low_memory=False)
#vegetables = ["cucumber", "tomato", "lettuce", "asparagus", # Replace this with the columns of dataframe
#            "potato", "wheat", "barley"]                   # except Gene IDs
#farmers = ["Farmer Joe", "Upland Bros.", "Smith Gardening", #Replace this with the Gene IDs
#           "Agrifun", "Organiculture", "BioGoods Ltd.", "Cornylee Corp."]
# Here we go...
# First let's lose the extraneous stuff (for viewer)
# like the contig and site
time_points=df.columns
genes=df["ID"]
#harvest = np.array([[0.8, 2.4, 2.5, 3.9, 0.0, 4.0, 0.0], # and this with the data
#                    [2.4, 0.0, 4.0, 1.0, 2.7, 0.0, 0.0],
#                    [1.1, 2.4, 0.8, 4.3, 1.9, 4.4, 0.0],
#                    [0.6, 0.0, 0.3, 0.0, 3.1, 0.0, 0.0],
#                    [0.7, 1.7, 0.6, 2.6, 2.2, 6.2, 0.0],
#                    [1.3, 1.2, 0.0, 0.0, 0.0, 3.2, 5.1],
#                    [0.1, 2.0, 0.0, 1.4, 0.0, 1.9, 6.3]])

data=df[df.columns.difference(["ID"])].to_numpy()

fig, ax = plt.subplots()
#im = ax.imshow(harvest)
im=ax.imshow(data, cmap=argstr.map,vmin=0.0,vmax=1.0)
# We want to show all ticks...
#ax.set_xticks(np.arange(len(farmers)))
ax.set_xticks(np.arange(len(time_points)-1))
#ax.set_yticks(np.arange(len(vegetables)))
ax.set_yticks(np.arange(len(genes)))
# ... and label them with the respective list entries
#ax.set_xticklabels(farmers)
ax.set_xticklabels(time_points[1:],fontsize=4)
#ax.set_yticklabels(vegetables)
ax.set_yticklabels(genes,fontsize=4)

# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

# Loop over data dimensions and create text annotations.
#for i in range(len(vegetables)):
#    for j in range(len(farmers)):
#        text = ax.text(j, i, harvest[i, j],
#                       ha="center", va="center", color="w")
for i in range(len(genes)):
    for j in range(len(time_points)-1):
        text = ax.text(j, i, data[i, j],
                       ha="center", va="center", color="w",fontsize=4)

#ax.set_title("Harvest of local farmers (in tons/year)")
ax.set_title("Methylation variation across time points")
fig.tight_layout()
plt.savefig(argstr.output,dpi=300,format="png")
