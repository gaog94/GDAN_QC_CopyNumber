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

### Significantly Altered Focal Peaks Analysis Code ###
def peakgene_overlaps(combos, same_genes, normalize=False):
    """Count the number of genes that overlap when examing the hg19 & hg38
    GISTIC runs' focal peaks."""
    venn_numbers, gsu, gsi = [], [], []
    for coh, ad in combos:
        print(coh)
        # put all significant genes in a list
        fnames = [hg19_dir+'/''+coh+ad+'genes.conf_99.txt', hg38_dir+'/'+coh+ad+'genes.txt']
        df38 = pd.read_csv(fnames[0], index_col=0).drop(['q value','residual q value','wide peak boundaries'], sep='\t')
        df19 = pd.read_csv(fnames[1], index_col=0).drop(['q value','residual q value','wide peak boundaries'], sep='\t')
        g_38 = set([x for col in df38.columns for x in df38[col].dropna()]) & same_genes
        g_19 = set([x for col in df19.columns for x in df19[col].dropna()]) & same_genes
        intersect, union = g_38 & g_19, g_38 | g_19
        gsu.append(union)
        gsi.append(intersect)
        if normalize:
            venn_numbers.append([len(g_19-intersect)/len(union),len(intersect)/len(union), len(g_38-intersect)/len(union)])
        else:
            venn_numbers.append([len(g_19-intersect),len(intersect), len(g_38-intersect)])

    index = [x[0]+'_'+x[1][1:-1] for x in combos]
    return pd.DataFrame(venn_numbers, index=index, columns=['hg19 only','Intersection','hg38 only'])

def plot_peakgene_overlaps(combos, same_genes, write=False):
    """Visualize the results of peakgene_overlaps function in bargraph form."""
    df_out = peakgene_overlaps(combos, same_genes, normalize=False)
    df_d, df_a = df_out[df_out.index.str.split('_').str[-1] == 'del'], \
                    df_out[df_out.index.str.split('_').str[-1] == 'amp']
    for x in zip((df_d, df_a), ('Deletion Peak Memberships', 'Amplification Peak Memberships')):
        x[0].index = x[0].index.str.split('_').str[0]
        x[0].plot.bar(stacked=True, color=['#af8dc3', '#f7f7f7', '#7fbf7b'], linewidth=1, edgecolor='k')
        plt.gca().set_xticklabels(x[0].index, rotation=90)
        plt.title(x[1], fontsize=18)
        plt.gcf().set_size_inches(10,8)
        sns.despine()
        plt.savefig(x[1].split(' ')[0]+'_peakMemberships.pdf', bbox_inches='tight')
        plt.savefig(x[1].split(' ')[0]+'_peakMemberships.png', bbox_inches='tight')
        plt.clf()
    if write:
        df_out.to_csv('VennStats_focalpeaks.tsv', sep='\t')


if __name__ == "__main__":

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
    ads = ['/amp_', '/del_']
    combos = [(c, a) for c in cohs for a in ads]

    # grab list of genes present in both hg19 & hg38 (use CHOL as a representative cancer type)
    df_hg38 = pd.read_csv(hg38_dir+'/CHOL/all_thresholded.by_genes.txt', index_col=0, usecols=[0,1], sep='\t')
    df_hg19 = pd.read_csv(hg19_dir+'/CHOL/all_thresholded.by_genes.txt', index_col=0, usecols=[0,1], sep='\t')
    same_genes = set(df_hg38.index) & set(df_hg19.index)

    plot_peakgene_overlaps(combos, same_genes, write=True)
