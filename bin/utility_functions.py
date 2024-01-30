#!/usr/bin/env python3

"""
Description: utility functions for plotting simple TCR repertoire statistics

Authors: Domenick Braccia
"""

## import packages
import time
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def TicTocGenerator():
    # Generator that returns time differences
    ti = 0           # initial time
    tf = time.time() # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti # returns the time difference

TicToc = TicTocGenerator() # create an instance of the TicTocGen generator

# This will be the main function through which we define both tic() and toc()
def toc(tempBool=True):
    # Prints the time difference yielded by generator instance TicToc
    tempTimeInterval = next(TicToc)
    if tempBool:
        print( "Elapsed time: %f seconds.\n" %tempTimeInterval )

def tic():
    # Records a time in TicToc, marks the beginning of a time interval
    toc(False)

# Defining sample comparison functions
def jaccard_index(sample1, sample2):
    set1 = set(sample1)
    set2 = set(sample2)
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union

def sorensen_index(sample1, sample2):
    set1 = set(sample1)
    set2 = set(sample2)
    intersection = len(set1.intersection(set2))
    return 2 * intersection / (len(set1) + len(set2))

def morisita_horn_index(dfs, sample1, sample2):
    # create sets of amino acid sequences
    set1 = set(dfs[sample1]['aminoAcid'])
    set2 = set(dfs[sample2]['aminoAcid'])
    # identify union of sets
    union = set1.union(set2)
    # loop through union of aa sequences and calculate morisita index between sample1 and sample2
    products=[]
    for aa in union:
        # get counts of aa sequences in sample1 and sample2
        if aa not in set(dfs[sample1]['aminoAcid']):
            n1i = 0
        else:
            n1i = dfs[sample1].where(dfs[sample1]['aminoAcid'] == aa)['read_count'].dropna().values[0]

        if aa not in set(dfs[sample2]['aminoAcid']):
            n2i = 0
        else:
            n2i = dfs[sample2].where(dfs[sample2]['aminoAcid'] == aa)['read_count'].dropna().values[0]
            
        product = n1i * n2i
        products.append(product)

    # calculate simpson index values for sample1 and sample2
    s1_si = sum([(count/sum(dfs[sample1]['read_count']))**2 for count in dfs[sample1]['read_count']])
    s2_si = sum([(count/sum(dfs[sample2]['read_count']))**2 for count in dfs[sample2]['read_count']])

    numerator = 2 * sum(products)
    denominator = (s1_si + s2_si) * (len(set1)*len(set2))
    return numerator / denominator

# def morisita_horn_index(sample1, sample2):
#     N1 = sum(sample1)
#     N2 = sum(sample2)
#     sum_n1i_n2i = sum([n1i * n2i for n1i, n2i in zip(sample1, sample2)])
#     sum_n1i_sq = sum([n1i**2 for n1i in sample1])
#     sum_n2i_sq = sum([n2i**2 for n2i in sample2])
#     return 2 * sum_n1i_n2i / ((sum_n1i_sq + sum_n2i_sq) * (N1 + N2))


#### =========================== LEGACY CODE ============================== ####

# def plot_timecourse2(df, x_col, y_col, patient_col):
#     # Create a list of colors for the scatter plot points
#     colors = []
#     for timepoint in df[x_col]:
#         if timepoint == 'Base':
#             colors.append('blue')
#         elif timepoint == 'Post':
#             colors.append('orange')

#     # Create a scatter plot of the data with the specified colors
#     plt.scatter(df[x_col], df[y_col], c=colors)

#     # Find the indices of the Base timepoints
#     base_indices = df[df[x_col] == 'Base'].index

#     # Iterate over the Base timepoints and plot lines to the corresponding Post timepoints
#     for base_idx in base_indices:
#         # Get the x and y coordinates of the Base timepoint
#         base_x, base_y = df.loc[base_idx, [x_col, y_col]]

#         # Find the index of the corresponding Post timepoint (if it exists)
#         post_idx = df[(df[patient_col] == df.loc[base_idx, patient_col]) & 
#                       (df[x_col] == 'Post')].index
#         if len(post_idx) > 0:
#             # Get the x and y coordinates of the Post timepoint
#             post_x, post_y = df.loc[post_idx[0], [x_col, y_col]]

#             # Plot a line between the Base and Post timepoints
#             plt.plot([base_x, post_x], [base_y, post_y], color='black')

#     # Add labels and title to the plot
#     plt.xlabel(x_col)
#     plt.ylabel(y_col)

# def plt_combined (df, x_col, y_col, patient_col='patient_id'):
#     fig, axs = plt.subplots(1, 2, figsize=(10, 5))
#     fig.tight_layout(pad=2)

#     ## box/strip plots overlaid
#     sns.boxplot(data=df, x=x_col, y=y_col, showfliers=False, ax=axs[0], color='white')
#     sns.stripplot(data=df, x=x_col, y=y_col, color='black', size=4, ax=axs[0])
#     axs[0].set(xlabel='', ylabel='')

#     ## line plot using plot_timecourse2()
#     plot_timecourse2(df, x_col=x_col, y_col=y_col, patient_col=patient_col)
#     axs[1].set(ylabel='', xlabel='')

#     ## figure adjustments
#     fig.subplots_adjust(top=0.90)

#     ## add title
#     titles = {
#         'num_clones': 'Number of unique T cell Clones',
#         'clonality': 'Clonality of T cell repertoire',
#         'simpson_index_corrected': 'Corrected Simpson index',
#         'pct_prod': 'Percentage of productive TCR sequences',
#         'cdr3_avg_len': 'Average Length of CDR3 sequences'
#     }
#     fig.suptitle(titles[y_col], fontsize=16)

#     ## save the plot
#     print('saving figure as: ' + y_col + '.png')
#     plt.savefig(y_col + '.png')

#     ## show the plot
#     plt.show()