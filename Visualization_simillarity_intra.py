#Short program to visualize the aa simillarity
#Author:Tarikul Islam Milon
#Created on: 29/03/2023
#Input files: 1) Simillarity_between_amino_acid_pairs_{aa}_{aa}.txt

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import itertools
import os
plt.figure(figsize=(12, 7.195), dpi=300)

# Enter the amino acids for visualizing the Box plot
#aa_list=['PTR', 'TPO', 'SEP']
#aa_list=['GLY','ALA','VAL','LEU','ILE','PRO','PHE','TYR','TRP','SER','THR','CYS','MET','HIS','LYS','ARG','ASP','GLU','ASN','GLN']
#aa_pair=itertools.combinations(aa_list,2)
aa_pair=['HIS','TRP']
labels=['Mirror','No_Mirror']
#for l in aa_list:
    #aa_pair.append((l,l)) #adding the aa itself pairs

   #Copy the folder path here
#aa_list=['VAL','TYR']
simillarity_dict={}
for aa in aa_pair:
    for label_ in labels:
        Folder_Path = f"/ddnA/work/wxx6941/TSR/code/code/Tarikul_code/Drug/{aa.lower()}/{label_}"

        f1 = open(f'{Folder_Path}/Simillarity_between_drug_pairs_{aa}_{aa}.txt',  'r')  # reads the .txt file for aa simillarity
        lines = f1.readlines()
        lines.remove(lines[0])
        Label=f'{aa}_{label_}'
        for line in lines:
            data = float(line.split()[2])
            if Label not in simillarity_dict:
                simillarity_dict[Label] = []
                simillarity_dict[Label].append(data)
            else:
                simillarity_dict[Label].append(data)

my_data=dict([(k,pd.Series(v)) for k,v in simillarity_dict.items()])

df0=pd.DataFrame(my_data)
df1=pd.melt(df0)


sns.set(font_scale=1.8)
#sns.set(rc={"figure.dpi":300, 'savefig.dpi':300})
#sns.set_theme(rc={'figure.dpi': 300, 'figure.figsize': (1, 1)}, font_scale=0.35)
sns.stripplot(x="variable", y="value", data=df1,jitter=1,hue='variable',palette='deep',legend=False, size=3.5,dodge=True)
sns.boxplot(x="variable", y="value", data=df1,showfliers = False,saturation=10,showmeans=True,width=0.5, linewidth = 2,meanprops={"marker":"s","markerfacecolor":"white", "markeredgecolor":"blue"})

plt.xticks(rotation=90)
#plt.title("Amino acids simillarity visualization")
plt.xlabel(None)
plt.ylabel('Similarity (%)',fontsize=25)


plt.savefig('Visualization_aa_similarity.png',dpi=300,bbox_inches='tight', facecolor='w', edgecolor='w', pad_inches=0.5)
#Display the plot
plt.show()
print('Completed successfully')



