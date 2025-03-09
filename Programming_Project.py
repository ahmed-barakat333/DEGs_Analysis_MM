'''
CIT 656 - Programming for Bioinformatics
Spring 2024
Final project
Idea Three
Students' names: Mohamed Ayman Mohamed, Khaled Mahmoud Saad,
Hossameldin Ibrahim Sharaf, Ahmed Barakat Ibrahim Soltan
'''

############################################################################################################
############################################################################################################


'''
The disease is multiple myeloma (a type of cancer that affects bone marrow plasma cells).

The dataset was downloaded from Gene Expression Omnibus (GEO) database, accession number GSE175384.
(URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175384)

Technology used is Bulk RNA-seq/ Illumina HiSeq 2000.

This dataset consists of 32 healthy samples and 32 diseased samples according to annotations provided by authors.

Raw dataset was processed where gene expression values were normalized and filtered using a custom
R script written as part of another project.
'''

## read in the data and perform statistical analysis

# take input from user to determine whether paired or independent t-test should be used
# check for valid input

while True:
    test_type = input('Which method should be used to analyse the data? \n'
                      'please choose between "paired" or "independent" or "wilcoxon"\n'
                      '')
    test_type = test_type.lower()
    if test_type not in ("independent", "paired", "wilcoxon"):
        print("Not an appropriate choice.\n")
    else:
        break


# read in dataset
import numpy as np
import pandas as pd
my_data = pd.read_csv('/Users/ahmedbarakat/Work/Diploma/Programming/Assignments/GSE175384_Filtered(ID,Counts).Normalized.Annotated.Data_voom.csv', sep=',', header=0)
            # replace with your filepath
print('\nHead of the dataset\n')
print(my_data.head()) # print head of the dataset

# rearrange dataframe columns based on sample annotations so that similar samples lie next to each other
cols= ['entrezgene_id','id','RNAseq33', 'RNAseq34', 'RNAseq35', 'RNAseq36', 'RNAseq37', 'RNAseq38', 'RNAseq39', 'RNAseq40', 'RNAseq41', 'RNAseq42',
           'RNAseq43', 'RNAseq44', 'RNAseq45', 'RNAseq46','RNAseq47', 'RNAseq48', 'RNAseq49', 'RNAseq50', 'RNAseq51', 'RNAseq52',
           'RNAseq53', 'RNAseq54', 'RNAseq55', 'RNAseq56', 'RNAseq57', 'RNAseq58', 'RNAseq59','RNAseq60', 'RNAseq61', 'RNAseq62',
           'RNAseq63', 'RNAseq64', # healthy samples

           'RNAseq100', 'RNAseq101','RNAseq102','RNAseq103','RNAseq104','RNAseq105',
           'RNAseq106', 'RNAseq107', 'RNAseq108', 'RNAseq109', 'RNAseq110', 'RNAseq111',
           'RNAseq80', 'RNAseq81', 'RNAseq82','RNAseq83', 'RNAseq84', 'RNAseq85', 'RNAseq86',
            'RNAseq87', 'RNAseq88', 'RNAseq89', 'RNAseq90', 'RNAseq91', 'RNAseq92', 'RNAseq93',
            'RNAseq94', 'RNAseq95', 'RNAseq96','RNAseq97', 'RNAseq98', 'RNAseq99' ] # diseased samples
my_data = my_data[cols]


# reshape dataframe from wide to long format to be suitable for statistical analysis
my_data_long = pd.melt(my_data, id_vars='id', value_vars=['RNAseq33', 'RNAseq34', 'RNAseq35', 'RNAseq36', 'RNAseq37', 'RNAseq38', 'RNAseq39', 'RNAseq40', 'RNAseq41', 'RNAseq42',
           'RNAseq43', 'RNAseq44', 'RNAseq45', 'RNAseq46','RNAseq47', 'RNAseq48', 'RNAseq49', 'RNAseq50', 'RNAseq51', 'RNAseq52',
           'RNAseq53', 'RNAseq54', 'RNAseq55', 'RNAseq56', 'RNAseq57', 'RNAseq58', 'RNAseq59','RNAseq60', 'RNAseq61', 'RNAseq62',
           'RNAseq63', 'RNAseq64', # healthy samples

           'RNAseq100', 'RNAseq101','RNAseq102','RNAseq103','RNAseq104','RNAseq105',
           'RNAseq106', 'RNAseq107', 'RNAseq108', 'RNAseq109', 'RNAseq110', 'RNAseq111',
           'RNAseq80', 'RNAseq81', 'RNAseq82','RNAseq83', 'RNAseq84', 'RNAseq85', 'RNAseq86',
            'RNAseq87', 'RNAseq88', 'RNAseq89', 'RNAseq90', 'RNAseq91', 'RNAseq92', 'RNAseq93',
            'RNAseq94', 'RNAseq95', 'RNAseq96','RNAseq97', 'RNAseq98', 'RNAseq99' ],
                       var_name='sample_name', value_name='value')


# define a function to fill in a new column to describe each sample group
healthy_samples = ['RNAseq33', 'RNAseq34', 'RNAseq35', 'RNAseq36', 'RNAseq37', 'RNAseq38', 'RNAseq39', 'RNAseq40', 'RNAseq41', 'RNAseq42',
                    'RNAseq43', 'RNAseq44', 'RNAseq45', 'RNAseq46','RNAseq47', 'RNAseq48', 'RNAseq49', 'RNAseq50', 'RNAseq51', 'RNAseq52',
                    'RNAseq53', 'RNAseq54', 'RNAseq55', 'RNAseq56', 'RNAseq57', 'RNAseq58', 'RNAseq59','RNAseq60', 'RNAseq61', 'RNAseq62',
                    'RNAseq63', 'RNAseq64']

def sample_group(sample_name):
    if any(x in sample_name for x in healthy_samples):
        return "Healthy"
    else:
        return "Diseased"

my_data_long['sample_group'] = my_data_long['sample_name'].map(sample_group)

print('\nHead of the long format dataset\n')
print(my_data_long.head())

# perform statistical analysis i.e., hypothesis testing
from scipy.stats import ttest_ind, ttest_rel, wilcoxon,false_discovery_control

# determine test type based on user input
if test_type == 'independent':
    res = my_data_long.groupby('id').apply(lambda x:
                                      ttest_ind(x[x['sample_group'] == "Healthy"]["value"],
                                                x[x['sample_group'] == "Diseased"]["value"]))
elif test_type == 'paired':
    res = my_data_long.groupby('id').apply(lambda x:
                                      ttest_rel(x[x['sample_group'] == "Healthy"]["value"],
                                                x[x['sample_group'] == "Diseased"]["value"]))

elif test_type == 'wilcoxon':
    res = my_data_long.groupby('id').apply(lambda x:
                                      wilcoxon(x[x['sample_group'] == "Healthy"]["value"],
                                               x[x['sample_group'] == "Diseased"]["value"]))

res = res.reset_index()

# extract output of statistical analysis
res[['t_value', 'p_value']] = pd.DataFrame(res[0].tolist(), index=res.index)

res = res[['id', 't_value', 'p_value']]

# multiple test correction of p_value
res['adj.p_value']= false_discovery_control(res['p_value'])

print('\nResults of hypothesis testing\n')
print(res)

# calculate log2 fold change
l2fc = pd.DataFrame(my_data_long.groupby('id').apply(lambda x:
                                       np.log2(np.nan_to_num(np.divide((x[x['sample_group']=="Healthy"]["value"]).mean(),(x[x['sample_group']=="Diseased"]["value"]).mean()), nan=1))))

print('\nResults of fold change\n')
print(l2fc.head())

# merge log2 fold change results with statistical analysis results
res_all = pd.merge(res, l2fc, on='id')
res_all.rename(columns={res_all.columns[-1]: 'l2fc'}, inplace=True)


# define a function to fill in a new column to describe gene expression change
def change_type(l2fc):
    if l2fc > 0:
        return "over-expressed"
    elif l2fc < 0:
        return "down-expressed"
    elif l2fc == 0:
        return "equally-expressed"

res_all['change_type'] = res_all['l2fc'].map(change_type)

############################################################################################################
############################################################################################################

## print and save results of differential gene expression

print('\nResults of hypothesis testing and fold change\n')
print(res_all.head())

# filter results based on adj.p_value or/and l2fc
res_all_pvalue = res_all[res_all['adj.p_value'] < 0.05]
res_all_l2fc = res_all[(res_all['l2fc'] > 0) | (res_all['l2fc'] < 0)]
res_all_pvalue_l2fc = res_all[(res_all['l2fc'] > 0) & (res_all['adj.p_value'] < 0.05)]._append(res_all[(res_all['l2fc'] < 0) & (res_all['adj.p_value'] < 0.05)])

print('\nThe number of differentially expressed genes (DEGs) based on adjusted p_value is %d' %(len(res_all_pvalue)))
print('\nHead of differentially expressed genes (DEGs) based on adjusted p_value is\n')
print(res_all_pvalue['id'].head())

print('\nThe number of differentially expressed genes (DEGs) based on log fold change is %d' %(len(res_all_l2fc)))
print('\nHead of differentially expressed genes (DEGs) based on log fold change is\n')
print(res_all_l2fc['id'].head())

print('\nThe number of differentially expressed genes (DEGs) based on adjusted p_value and log fold change is %d' %(len(res_all_pvalue_l2fc)))
print('\nHead of differentially expressed genes (DEGs) based on adjusted p_value and log fold change is\n')
print(res_all_pvalue_l2fc['id'].head())

# print top 10 differentially expressed genes
top_ids = res_all_pvalue.sort_values(by=["adj.p_value","l2fc"], ascending=[False,False])
top_ids = top_ids.head(10)
print('\nTop ten differentially expressed genes (DEGs) are\n %s'%(top_ids['id'].to_list()))

# save results into file
res_all_pvalue_file = ('Hypothesis testing results' + "\n" + "\n" + res_all_pvalue.to_string(index=False))
res_all_l2fc_file = ('Fold change results' + "\n" + "\n" + res_all_l2fc.to_string(index=False))
res_all_pvalue_l2fc_file = ('Hypothesis testing and fold change results' + "\n" + "\n" + res_all_pvalue_l2fc.to_string(index=False))

destination_file = open('/Users/ahmedbarakat/Work/Diploma/Programming/Assignments/final_project_degs.txt', 'w')
destination_file.write(res_all_pvalue_file + '\n' + '\n' + res_all_l2fc_file)


############################################################################################################
############################################################################################################

## plot results of differential gene expression (DEG)

# prepare the two sets of hypothesis testing and fold change analysis results
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
set1 = set(res_all_pvalue['id'])
set2 = set(res_all_l2fc['id'])

# plot  venn diagram
venn = venn2([set1, set2], ('Hypothesis testing', 'Fold change'))
plt.title('Venn diagram of hypothesis testing and fold change analysis results')

plt.show()

# volcano plot
# reference tutorial URL: https://ai.plainenglish.io/volcano-plots-9721d88fc3ef
import matplotlib as mpl

alpha = 0.05
plt.figure(figsize=(8, 6))
plt.scatter(res_all['l2fc'], -np.log10(res_all['adj.p_value']), color='blue', alpha=0.5)
plt.axhline(-np.log10(alpha), color='red', linestyle='--', linewidth=1)
plt.xlabel('Log2 Fold Change')
plt.ylabel('-log10(p-value)')
plt.title('Volcano plot')
plt.grid(True)
plt.show()

# heatmap plot of gene expression values of top 66 genes according to adj.p_value
# reference tutorial URL: https://medium.com/@SeanAT19/how-to-create-heatmaps-using-matplotlib-pyplot-db9ff94d25e8

# get DEGs with highest log2 fold change e.g., 66 to match no. of samples
heatmap = my_data[my_data['id'].isin(res_all.nlargest(66, ['l2fc'])['id'])]
fig, ax = plt.subplots()

# define colors and bounds of heatmap
colors = ['blue', 'white', 'red']
bounds = [0,1,2,3]

cmap = mpl.colors.ListedColormap(colors)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
im = ax.imshow(heatmap.iloc[1:, 2:], interpolation='none', cmap=cmap, norm=norm)
plt.colorbar(im)

ax.set_yticks(np.arange(len(heatmap['id'])), labels=heatmap['id'])
ax.set_xticks(np.arange(len(list(heatmap.columns[2:]))), labels=list(heatmap.columns[2:]))
plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
         rotation_mode="anchor")
plt.xlabel('Sample id')
plt.ylabel('Gene id')
plt.title('Heatmap of top differentially expressed genes (DEGs)')
plt.show()

############################################################################################################
############################################################################################################

## construct a protein-protein interaction network using DEGs
# get PPI data
# URL: https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/
interaction_data = pd.read_csv('/Users/ahmedbarakat/Work/Diploma/Intro/Assignments/hippie_current.txt', sep="\t", header=None, on_bad_lines='skip',
                               usecols=[0, 2, 4], low_memory=True)
                # replace with your filepath

# extract certain columns
interaction_data.rename(columns={0: 'head', 2: 'tail', 4: 'weight'}, inplace=True)
interaction_data['head'] = interaction_data['head'].str.split('_').str[0]
interaction_data['tail'] = interaction_data['tail'].str.split('_').str[0]

# filter PPI dataset by edge weight >= 0.9
interaction_data = interaction_data.loc[interaction_data['weight'] >= 0.9]

# filter PPI dataset by DEG list
interaction_data_filtered = interaction_data[interaction_data['head'].isin(res_all_pvalue['id'].tolist())]

# remove self interacting proteins
interaction_data_filtered = interaction_data_filtered[interaction_data_filtered['head'] != interaction_data_filtered['tail']]

# construct the network
import networkx as nx
G=nx.from_pandas_edgelist(interaction_data_filtered, "head", "tail")

# remove isolated nodes and keep the biggest connected subnetwork
G.remove_nodes_from(list(nx.isolates(G)))
Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
G0 = G.subgraph(Gcc[0])

# plot the network
plt.figure(figsize=(40, 40))
nx.draw(G0, with_labels= True)
plt.savefig("network.png")

# analyze the network to get nodes with highest degree i.e., hub nodes
degrees = dict(G.degree())

# define a function to get number of nodes with highest degree
def top_x_nodes(degrees, x):
    sorted_nodes = sorted(degrees, key=degrees.get, reverse=True)
    top_nodes = sorted_nodes[:x]
    return top_nodes

# choose the number of nodes to be returned
x = 15
mm_network = top_x_nodes(degrees, x)

# print top 10 nodes and their degrees
print(f"\nTop {x} nodes with highest degree:\n")
for node in mm_network:
    print(f"Node {node}: Degree = {degrees[node]}")