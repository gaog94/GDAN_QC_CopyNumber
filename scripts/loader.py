#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 17:28:54 2018

@author: galengao
"""
import pandas as pd

### Helper Function to Load in the Data ###
def load_data(hg38_dir, hg19_dir, coh, thresh=False):
    """Load in the hg38 and hg19 gistic thresholded data. Assume GISTIC runs
    for each tumor type live in a parent directory (hg38_gistic or hg19_gistic)
    one level up from this script."""
    if thresh:
        hg38 = hg38_dir+'/'+coh+'/all_thresholded.by_genes.txt'
        hg19 = hg19_dir+'/'+coh+'/all_thresholded.by_genes.txt'
        hg38drops = ['Cytoband', 'Locus ID']
    else:
        hg38 = hg38_dir+'/'+coh+'/all_data_by_genes.txt'
        hg19 = hg19_dir+'/'+coh+'/all_data_by_genes.txt'
        hg38drops = ['Cytoband', 'Gene ID']

    df_hg19 = pd.read_csv(hg19, index_col=[0]).drop(['Cytoband', 'Locus ID'], axis=1, sep='\t')
    df_hg38 = pd.read_csv(hg38, index_col=[0]).drop(hg38drops, axis=1, sep='\t')

    same_samps = list(set(df_hg38.columns) & set(df_hg19.columns))
    same_genes = list(set(df_hg38.index) & set(df_hg19.index))
    print(coh, len(same_genes), len(same_samps))
    return df_hg38[same_samps].T[same_genes], df_hg19[same_samps].T[same_genes]

    return df_hg38, df_hg19
