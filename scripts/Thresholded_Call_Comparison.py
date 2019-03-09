#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 17:28:54 2018

@author: galengao
"""
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import loader

### Thresholded Copy Number Values Analysis Code ###
def thresholded_value_comparison(df_hg38, df_hg19, metric='hamming'):
    """Compare -2,-1,0,1,2 gene-level thresholded calls. metric can be either
    hamming (number of discrepancies in each gene) or manhattan (sum of
    'distances' between each gene so a 1 to -1 change is 2). Returns a vector
    of each gene's metric."""
    out = []
    for i, g in enumerate(df_hg38.columns):
        if metric == 'hamming':
            out.append(sum(df_hg19[g] != df_hg38[g])/len(df_hg19))
        elif metric == 'manhattan':
            out.append(sum(abs((df_hg19[g] - df_hg38[g]))))
    return pd.DataFrame(out, index=df_hg38.columns)

def sequential_cohort_test_thresholded_values(cohs):
    """Compare thresholded gene-level calls for input tumor types."""
    df_out = pd.DataFrame([])
    for coh in cohs:
        df_hg38, df_hg19 = loader.load_data(hg38_dir, hg19_dir, coh, thresh=True)
        df_results = thresholded_value_comparison(df_hg38, df_hg19, metric='hamming')
        df_results.columns = [coh]
        df_out = df_out.join(df_results, how='outer')

    df_out.to_csv('../readout/DiscordantSampleFractions_perGene_perCohort_thresholdedCalls.tsv', sep='\t')
    return df_out

def plot_fractionDisagreements_perCohort(cohs):
    """Visualize fraction of samples with disagreements in thresholded copy
    number for each gene. Run sequential_cohort_test_thresholded_values()
    before this function."""
    # Read in data written by sequential_cohort_test_thresholded_values
    df = sequential_cohort_test_thresholded_values(cohs)
    df_box = pd.melt(df.reset_index(), id_vars='Gene Symbol').set_index('Gene Symbol')
    df_box.columns = ['Tumor Type', 'Fraction of Samples with Disagreements']
    dft = df.T
    dft['med_degenerates'] = df.median(axis=0)
    boxorder = dft.sort_values('med_degenerates', axis=0).index

    # read in copy number burden data (requires aneuploidy RecurrentSCNA calls)
    df_cn = pd.read_csv(aneuploidy_call_file, index_col=0, usecols=[0,1,2,16], sep='\t')
    coh_medians = [int(np.median(df_cn[df_cn['Type']==x]['RecurrentSCNA'].dropna())) for x in df_cn.Type.unique()]
    df_med = pd.DataFrame(coh_medians, index=df_cn.Type.unique(), columns=['med'])

    # plot it out
    pal = sns.color_palette('Blues', max(df_med.med)-min(df_med.med)+1)
    my_pal = {c: pal[df_med.at[c,'med']] for c in df_med.index}
    g = sns.boxplot(x=df_box.columns[0], y=df_box.columns[1], data=df_box, \
                    order=boxorder, fliersize=1, palette=my_pal, linewidth=0.5)
    newxticks = [x+' ('+str(df_med.loc[x]['med'])+')' for x in boxorder]
    g.set_xticklabels(newxticks, rotation=90)
    plt.ylabel('Fraction with Disagreements', fontsize=12)
    sns.despine()
    plt.gcf().set_size_inches((8,3))
    plt.savefig('2_thresholdedCN_boxplot.pdf', bbox_inches='tight')
    plt.savefig('2_thresholdedCN_boxplot.png', bbox_inches='tight')


if __name__ == "__main__":
    # parse cmd-line args
    if len(sys.argv) == 4:
      hg38_dir = sys.argv[1]
      hg19_dir = sys.argv[2]
      aneuploidy_call_file = sys.argv[3]
    else:
      hg38_dir = '../hg38_gistic/'
      hg19_dir = '../hg19_gistic/'
      aneuploidy_call_file = '../PANCAN_armonly_ASandpuritycalls_092817_xcellcalls.txt'
      print('No GISTIC output directories provided. Assuming they exist at ../hg19_gistic/ and ../hg38_gistic/')
      print('PanCancer Aneuploidy AWG list of number of focal alterations also not provided. Assuming it exists at ../PANCAN_armonly_ASandpuritycalls_092817_xcellcalls.txt')

    # set up the tumor types we want to analyze
    cohs = ['ACC','BLCA','CESC','CHOL','COAD','DLBC','ESCA','GBM', 'HNSC','KICH',\
            'KIRC','KIRP','LAML','LGG','LIHC','LUAD','LUSC','OV','PAAD','PCPG',\
            'PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM']

    plot_fractionDisagreements_perCohort(cohs)
