#!/usr/bin/env python
# coding: utf-8

# ## Generating intensity table for DHS sites
# 
# - Give an overview of data samples, with metadata
# - Run coverage_table in Gimme Motif
# 
# #### Loading the important dependencies

# In[1]:


get_ipython().run_line_magic('load_ext', 'nb_black')
get_ipython().run_line_magic('reload_ext', 'nb_black')


# In[2]:


import os
import pandas as pd
import subprocess as sp
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import seaborn as sns
import warnings
import qnorm

from sklearn.preprocessing import StandardScaler  # For scaling dataset
from sklearn.cluster import (
    KMeans,
    AgglomerativeClustering,
    AffinityPropagation,
)  # For clustering
from sklearn.decomposition import PCA


# #### Setting file paths

# In[3]:


## General locations
exp_path = "/scratch/snabel/GEODatasets_CMenhancers"
bam_path = "/scratch/snabel/GEODatasets_CMenhancers/results/testset/"

samplesheet = "/scratch/snabel/GEODatasets_CMenhancers/samplesheet_all.tsv"
metadata = "/scratch/snabel/GEODatasets_CMenhancers/metadata.tsv"

## Files for coverage_table
dhs_file = "/scratch/snabel/GEODatasets_CMenhancers/Meuleman_DHSs_2019/ENCFF503GCK_summit_10k.bed"
dhs_set = "10k"
coveragetable = "/scratch/snabel/GEODatasets_CMenhancers/10k_2000bp_window_coverage_table_logt.txt"
rerun_coverage_table = False


# #### Setting parameters

# In[4]:


window = "2000"
marks = ["H3K27ac", "H3K4me3"]
colors = {"H3K27ac": "#1F78B4", "H3K4me3": "#33A02C"}


# ### Selecting the metadata for the samples of interest
# Multiple samples did not run for the SRA dump, did not have an SRA entry or are an input sample (only needed earlier on in the workflow), therefore from the large metadata, only the samples used for analysis will be selected.

# In[5]:


samplesh_df = pd.read_table(samplesheet, na_values=["NaN"], keep_default_na=False)
# only keeping the files without an entry in the problems column
samplesh_df = samplesh_df[samplesh_df["problems"] == ""]
sample_ids = samplesh_df[["sample", "descriptive_name"]]

# large metatable with all samples + input
metadf = pd.read_table(metadata, na_values=["NaN"], keep_default_na=False)

# keep only the samples that were aligned
metadf_merge = pd.merge(
    left=sample_ids, right=metadf, on="sample", how="left", validate="one_to_one"
)
metadf_merge[:5]


# #### How many samples per source tissue are present for the marks of interest?

# In[6]:


metadf_merge.source_tissue.value_counts()
metadf_merge.target.value_counts()
crosstab_df = pd.crosstab(
    index=metadf_merge["source_tissue"], columns=metadf_merge["target"], margins=True
)
# Show the tissues of which the marks of interest are present
crosstab_df[marks]


# ## Running coverage_table
# Coverage_table needs as input the bam files and a bed file, the bed file defines the regions in which the reads will be counted per sample. In this case the DNase I hypersensibility sites (Meuleman et al. 2020) will be used as open regions to query. Selecting the bam files to run.

# In[14]:


bamfiles_path = [
    os.path.abspath(os.path.join(bam_path, p))
    for p in os.listdir(bam_path)
    if p.endswith((".bam"))
]

# generating the same list without file path
bamfiles = [file for file in os.listdir(bam_path) if file.endswith(".bam")]
bamfiles = [re.sub(".samtools-coordinate.bam", "", x) for x in bamfiles]
bamfiles = [re.sub("^GRCh38-", "", x) for x in bamfiles]


# Select from the metadata the samples used for running coverage_table.

# In[10]:


samples_meta = metadf_merge[np.isin(metadf_merge["descriptive_name"], bamfiles)]


# #### Running the coverage_table command in bash
# Coverage_table will also perform a log2 transformation.

# In[11]:


if rerun_coverage_table == True:
    print("running coverage table ")
    input_file = " ".join(bamfiles_path)
    filename = dhs_set + "_" + window + "bp_window_coverage_table_logt"
    output_file = filename + ".txt"
    log_file = filename + ".log"
    print(
        "using "
        + str(len(bamfiles_path))
        + " supplied files (check "
        + output_file
        + ")"
    )
    sp.check_call(
        f"coverage_table "
        f"-p {dhs_file} "
        f"-d {input_file} "
        f"> {exp_path}/{output_file} "
        f"-w {window} -l "
        f"2> {exp_path}/{log_file}",
        shell=True,
    )
else:
    output_file = coveragetable


# #### Load the coverage table per mark into a dataframe

# Generate clean sample names for the coverage_table columns using the bamfilenames, stripped of assembly name and filename extension. 

# In[15]:


# The first column of the coverage_table output, is always the location/region
cov_column_names = bamfiles
cov_column_names.insert(0, "region")

## Create a df with the loaded output file, containing 9 rows of information regarding the run
coverage_table_df = pd.read_table(
    f"{output_file}", skiprows=9, index_col=0, names=cov_column_names
)
coverage_table_df[:5]


# In[18]:


# Make a dictionary with a df per mark.
dict_of_dfs = {}

for mark in marks:
    mark_samples = samples_meta.loc[samples_meta["target"] == mark].descriptive_name
    dict_of_dfs[mark] = coverage_table_df[mark_samples]


# ### Generate a distribution plot per sample. 
# Before and after quantile normalization.

# In[20]:


warnings.filterwarnings('ignore')

plt.style.use("seaborn")
sns.set(rc={"figure.figsize": (23, 16)})
# how to make this grid variable in size, depending on the 
# amount of samples I have. (now only 15 samples, but I want to run for >100)
fig, axes = plt.subplots(3, 5)

plt_x_ax = 0
plt_y_ax = 0

for mark in marks: 
    samples = list(dict_of_dfs[mark])
    for sample in samples:
        sns.distplot(dict_of_dfs[mark].get(sample),
                 hist=True,
                 rug=False,
                 label=sample,
                 axlabel=(sample),
                 norm_hist=False,
                 color=colors[mark],
                 ax=axes[plt_y_ax,plt_x_ax],
        )
        if plt_x_ax < 4:
            plt_x_ax = plt_x_ax + 1
        elif plt_x_ax == 4:
            plt_x_ax = 0
            plt_y_ax = plt_y_ax + 1

    fig.savefig(f"{exp_path}/Distribution_graphs_notNorm.pdf")


# #### Quantile normalization on samples per mark

# In[21]:


dict_norm_dfs = {}
for mark in marks:
    df = dict_of_dfs[mark]
    norm_df = qnorm.quantile_normalize(df, axis=1, ncpus=20)
    dict_norm_dfs[mark] = norm_df


# In[22]:


warnings.filterwarnings("ignore")

plt.style.use("seaborn")
sns.set(rc={"figure.figsize": (23, 16)})
fig, axes = plt.subplots(3, 5)

plt_x_ax = 0
plt_y_ax = 0

for mark in marks:
    samples = list(dict_norm_dfs[mark])
    for sample in samples:
        sns.distplot(
            dict_norm_dfs[mark].get(sample),
            hist=True,
            rug=False,
            label=sample,
            axlabel=(sample),
            norm_hist=False,
            color=colors[mark],
            ax=axes[plt_y_ax, plt_x_ax],
        )
        if plt_x_ax < 4:
            plt_x_ax = plt_x_ax + 1
        elif plt_x_ax == 4:
            plt_x_ax = 0
            plt_y_ax = plt_y_ax + 1

    fig.savefig(f"{exp_path}/Distribution_graphs_qNorm.pdf")


# ### Exploratory visualization of the samples
# Using PCA and heatmap representations.

# In[24]:


# Generating a large df with both marks combined (which were normalized seperately)
merged_marks_df = pd.merge(
    left=dict_norm_dfs["H3K27ac"],
    right=dict_norm_dfs["H3K4me3"],
    on="region",
    how="inner",
    validate="one_to_one",
)


# #### Running PCA on the merged dataset

# In[25]:


# Transpose df for PCA; making columns features/attributes and rows samples
merged_marks_tdf = merged_marks_df.T
pca = PCA(n_components=2)
pcs = pca.fit_transform(merged_marks_tdf)
pc_df = pd.DataFrame(data=pcs, columns=["PC1", "PC2"])
pc_df = pd.concat([pd.DataFrame(merged_marks_df.columns), pc_df], axis=1)
pc_meta = samples_meta[["descriptive_name", "target", "source_tissue"]]
pc_df = pd.merge(
    left=pc_df,
    right=pc_meta,
    left_on=0,
    right_on="descriptive_name",
    how="left",
)
pc_df = pc_df.drop(columns=0)
pc_df[:5]


# In[26]:


print(
    "Explained variation per principal component: {}".format(
        pca.explained_variance_ratio_
    )
)


# In[27]:


# Plot PCA with source_tissue labelling
plt.figure()
plt.figure(figsize=(10, 10))
plt.xticks(fontsize=12)
plt.yticks(fontsize=14)
plt.xlabel("Principal Component - 1", fontsize=20)
plt.ylabel("Principal Component - 2", fontsize=20)
plt.title("Principal Component Analysis", fontsize=20)

source_tissues = samples_meta["source_tissue"].unique()
# still need to make the colorscheme of a variable length,
# depending on the amount of tissues.
# colors = cm.get_cmap("viridis", len(source_tissues)).colors
colors = ["r", "g", "b", "purple"]
for target, color in zip(source_tissues, colors):
    indicesToKeep = pc_df["source_tissue"] == target
    plt.scatter(
        pc_df.loc[indicesToKeep, "PC1"],
        pc_df.loc[indicesToKeep, "PC2"],
        c=color,
        s=50,
    )

plt.legend(source_tissues, prop={"size": 15})
plt.savefig(f"{exp_path}/PCA_PC1-2_tissuesource.pdf")


# In[28]:


# Plot PCA with mark labelling
plt.figure()
plt.figure(figsize=(10, 10))
plt.xticks(fontsize=12)
plt.yticks(fontsize=14)
plt.xlabel("Principal Component - 1", fontsize=20)
plt.ylabel("Principal Component - 2", fontsize=20)
plt.title("Principal Component Analysis", fontsize=20)

source_tissues = samples_meta["source_tissue"].unique()
colors = ["b", "g"]
for target, color in zip(marks, colors):
    indicesToKeep = pc_df["target"] == target
    plt.scatter(
        pc_df.loc[indicesToKeep, "PC1"],
        pc_df.loc[indicesToKeep, "PC2"],
        c=color,
        s=50,
    )

plt.legend(marks, prop={"size": 15})
plt.savefig(f"{exp_path}/PCA_PC1-2_marks.pdf")


# In[29]:


sns.clustermap(merged_marks_df)


# In[30]:


get_ipython().system('/vol/mbconda/snabel/anaconda3/condabin/conda list > {exp_path}/conda_env.txt')


# In[ ]:


from IPython.display import Javascript
display(Javascript('IPython.notebook.execute_cells_above()'))


# In[ ]:


## Use this example to make a PCA with size of sample as visual variable: 
# https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/scatter_demo2.html#sphx-glr-gallery-lines-bars-and-markers-scatter-demo2-py

