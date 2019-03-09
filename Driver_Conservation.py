#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 17:28:54 2018

@author: galengao
"""
import sys

import pandas as pd

### Conservation of Significant Copy Number Driver Events Analysis Code ###
def documented_driver_differences(fname):
    """Scan and analyze manually currated DocumentedDriverDifferences.txt file.
    Returns: 1) Number of driver genes called in both hg19 & hg38 GISTIC peaks
    2) Number of drivers missing in hg38 peaks that appeared in hg19 peaks  and
    3) Number of drivers present in hg38 peaks but absent from hg19 peaks."""
    # read in table of documented driver differences
    # (this table needs a manual curation to be generated)
    df = pd.read_table(fname, index_col=0)
    # process entries to have just yes/no calls (without parens & brackets)
    df['hg19?'] = df['present in hg19?'].str.strip(')').str.strip('(').str.strip('[').str.strip(']')
    df['hg38?'] = df['present in hg38?'].str.strip(')').str.strip('(').str.strip('[').str.strip(']')

    # number of documented drivers that match in hg19 & hg38
    matches = sum(df['hg19?'] == df['hg38?'])
    # number of documented drivers that are in hg19 but not hg38 & vice versa
    lostdrivers = len(df[(df['hg19?'] == 'yes') & (df['hg38?'] == 'no')])
    recovereddrivers = len(df[(df['hg19?'] == 'no') & (df['hg38?'] == 'yes')])

    # Return in order
    return matches, lostdrivers, recovereddrivers


if __name__ == "__main__":
    if len(sys.argv) == 0:
      fname = '../DocumentedDriverDifferences.txt'
      print('No file of documented driver differences provided. Assuming it exists at ../DocumentedDriverDifferences.txt')
    else:
      fname = sys.argv[1]

    print(documented_driver_differences(fname))
