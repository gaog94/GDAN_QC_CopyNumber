#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 17:28:54 2018

@author: galengao
"""
import numpy as np
import pandas as pd

from collections import Counter

import matplotlib.pyplot as plt
import seaborn as sns

import loader

### Raw Copy Number Values Analysis Code ###
def raw_value_comparison(coh, plot=False):
    """Return the average differences in raw copy number values between the
    gene-level calls in hg19 and hg38 for each gene for a given tumor type
    'coh.' If plot=True, plot the genes' differences in a histogram."""

    # load in the data
    df_38, df_19 = loader.load_data(hg38_dir, hg19_dir, coh, thresh=False)

    # compute average sample-by-sample differences for each gene
    df_s = df_38 - df_19
    avg_diff = {g:np.average(df_s[g]) for g in df_s.columns.get_level_values('Gene Symbol')}

    # take note of which genes are altered more than our threshold of 4*std
    results = []
    std = np.std([avg_diff[x] for x in avg_diff])
    for g in avg_diff:
        if avg_diff[g] > 4 * std:
            results.append([coh, 'Pos', g, avg_diff[g]])
        elif avg_diff[g] < -4 * std:
            results.append([coh, 'Neg', g, avg_diff[g]])

    if plot:
        plt.hist([avg_diff[x] for x in avg_diff], bins=1000)
        plt.title(coh, fontsize=16)
        plt.xlabel('Average CN Difference Between Alignments', fontsize=14)
        plt.ylabel('Genes', fontsize=14)
        sns.despine()
        plt.savefig('./genehists/'+coh+'_genehist.pdf')
        plt.savefig('./genehists/'+coh+'_genehist.png')
        plt.clf()

    return results

def sequential_cohort_test_raw_values(cohs, plot=False):
    """Sequentially compare raw gene-level calls for the given tumor types."""
    c_results = []
    for coh in cohs: # perform raw value comparison for each cohort
        c_results += raw_value_comparison(coh, plot=plot)

    # compile results together
    df_r = pd.DataFrame(c_results, columns=['Cohort', 'Direction', 'Gene', 'Difference'])
    gcount = Counter(df_r['Gene'])
    pos_gcount = Counter(df_r[df_r['Direction']=='Pos']['Gene'])
    neg_gcount = Counter(df_r[df_r['Direction']=='Neg']['Gene'])
    df = pd.DataFrame([gcount[x] for x in gcount], index=gcount.keys(), columns=['Count'])
    df['Count_pos'] = [pos_gcount[x] if x in pos_gcount else 0 for x in gcount]
    df['Count_neg'] = [neg_gcount[x] if x in neg_gcount else 0 for x in gcount]

    if plot: # write output
        plt.plot(np.sort([gcount[x] for x in gcount])[::-1], 'b-')
        plt.xlabel('Gene by Rank', fontsize=16)
        plt.ylabel('Number of Occurences', fontsize=16)
        sns.despine()
        plt.savefig('GeneDevianceDropoff.pdf')
        plt.savefig('GeneDevianceDropoff.png')
        df_r.to_csv('./genehists/LargestDifferences.tsv', sep='\t', index=False)
        df.to_csv('./genehists/LargestDifferenceGenes_ByCount.tsv', sep='\t', index=True)

if __name__ == '__main__':

    # parse cmd-line args
    if len(sys.argv) == 3:
      hg38_dir = sys.argv[1]
      hg19_dir = sys.argv[2]
    else:
      hg38_dir = '../hg38_gistic/'
      hg19_dir = '../hg19_gistic/'
      print('No GISTIC output directories provided. Assuming they exist at ../hg19_gistic/ and ../hg38_gistic/')

    # set up the tumor types we want to analyze
    cohs = ['ACC','BLCA','CESC','CHOL','COAD','DLBC','ESCA','GBM', 'HNSC','KICH',\
            'KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','PAAD','PCPG',\
            'PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM']

    sequential_cohort_test_raw_values(cohs, plot=True)
