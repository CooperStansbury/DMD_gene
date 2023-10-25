import pandas as pd
import numpy as np
import os
import sys

tmap = {
    'S1a' : 0,
    'S1b' : 0,
    'S2a' : 0,
    'S2b' : 0,
    'S3a' : 1,
    'S3b' : 1,
    'S4a' : 2,
    'S4b' : 2,
    'S5a' : 3,
    'S5b' : 3,
    'S6a' : 4,
    'S6b' : 4,
    'S7a' : 5,
    'S7b' : 5,
    'S8a' : 6,
    'S8b' : 6,
    'S9a' : 7,
    'S9b' : 7,
}

rmap = {
    'S1a' : "r1c",
    'S1b' : "r2c",
    'S2a' : "r1",
    'S2b' : "r2",
    'S3a' : "r1",
    'S3b' : "r2",
    'S4a' : "r1",
    'S4b' : "r2",
    'S5a' : "r1",
    'S5b' : "r2",
    'S6a' : "r1",
    'S6b' : "r2",
    'S7a' : "r1",
    'S7b' : "r2",
    'S8a' : "r1",
    'S8b' : "r2",
    'S9a' : "r1",
    'S9b' : "r2",
}

cmap = {
    'S1a' : "control",
    'S1b' : "control",
    'S2a' : "timecourse",
    'S2b' : "timecourse",
    'S3a' : "timecourse",
    'S3b' : "timecourse",
    'S4a' : "timecourse",
    'S4b' : "timecourse",
    'S5a' : "timecourse",
    'S5b' : "timecourse",
    'S6a' : "timecourse",
    'S6b' : "timecourse",
    'S7a' : "timecourse",
    'S7b' : "timecourse",
    'S8a' : "timecourse",
    'S8b' : "timecourse",
    'S9a' : "timecourse",
    'S9b' : "timecourse",
}

def getGeneLengths(gene_table_path, gene_names):
    """A function to get gene lengths from a gene_table file

    params:
        : gene_table_path (str): path to the geneTable.csv file(pipeline output)
        : gene_names (list of str): valid gene names
    """
    gf = pd.read_csv(gene_table_path)
    gf = gf[gf['gene_name'].isin(gene_names)]
    gf = gf[gf['Feature'] == 'gene']
    gf = gf[gf['gene_biotype'] == 'protein_coding'].reset_index(drop=False)
    gf = gf[['gene_name', 'Start', 'End']].drop_duplicates()
    gf['Length'] = gf['End'] - gf['Start']
    gf = gf.groupby(['gene_name'])['Length'].max().reset_index(drop=False)
    return gf


def TPM(df, gf, target=1e6, p=1000):
    """A function to compute TPM for each column if a dataframe 
    
    params:
        : df (pd.DataFrame): the data 
        : gf (pd.DataFrame): gene lengths
        : target (float): the normalized sum of reads
        : p (int): the gene length normalization factor, default is kilo
    """
    tpm = df.copy()
    for c in tpm.columns:
        reads_per_gene  = df[c] / (gf['Length'].to_numpy() / p)
        sequence_depth = df[c].sum() / target
        tpm[c] = reads_per_gene / sequence_depth
    return tpm


def CPM(df, target=1e6):
    """A function to compute CPM for each column if a dataframe 
    
    params:
        : df (pd.DataFrame): the data 
        : target (float): the normalized sum of reads
    """
    cpm = df.copy()
    for c in cpm.columns:
        cpm[c] = (df[c] / df[c].sum()) *  target
    return cpm