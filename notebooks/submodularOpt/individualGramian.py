# This function computes the Gramian if only individual sensors are selected and writes each one to a file.

# Imports
import pandas as pd
import numpy as np
from copy import deepcopy
import os
import sys
from importlib import reload
from scipy.stats import zscore
import scipy.io
from scipy import sparse
import leidenalg
import scipy
import textwrap

from pydmd import DMD

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# local imports
sys.path.append("../../python/")
sys.path.append("../")

import nb_util as nb
import utils as ut
reload(ut)

import hasnain_DMD
reload(hasnain_DMD)

def getC(n, idxs):
    C = np.zeros((len(idxs), n))
    for i in range(len(idxs)):
        C[i, int(idxs[i])] = 1
    return C

# Load data (taken exactly from Cooper)
data_path = f"/nfs/turbo/umms-indikar/shared/projects/cell_cycle/data/RNA_pipeline_ouputs/countMatrix/counts.raw.txt"
gene_path = f"/nfs/turbo/umms-indikar/shared/projects/cell_cycle/data/RNA_pipeline_ouputs/references/geneTable.csv"

""" Load the raw expression """
df = pd.read_csv(data_path, index_col=0)

# remove MT and ribosomal genes
all_genes = df.index.to_list()
mt_genes = [x for x in all_genes if x.startswith('MT-')]
rp_genes = [x for x in all_genes if x.startswith('RP')]

print(f"{df.shape=}")
df = df.drop(mt_genes) # drop MT genes
df = df.drop(rp_genes) # drop ribosomal genes
print(f"{df.shape=}")

# rewrite the list without MT genes
gene_names = df.index.to_list()

print(f"{len(all_genes)=} {len(mt_genes)=} {len(gene_names)=}")

""" Load gene lengths """
gf = nb.getGeneLengths(gene_path, gene_names)
print(f"{gf.shape=}")

target = 1e6
threshold = 0.5
rank = 7
tpm = nb.TPM(df, gf, target=target)

# get highly expressed genes
tpm_dist = tpm.mean(axis=1)
mask = (tpm_dist > threshold)
high_exp_genes = tpm_dist[mask].index.to_list()    

# filter and convert to fold changes
d = tpm[tpm.index.isin(high_exp_genes)]
dmd_data = nb.data2DMD(d) 
print(f"{dmd_data.shape=}")

dmd_res = hasnain_DMD.dmd(dmd_data, rank=rank)

print(f"{dmd_res['A'].shape=}")
print(f"{dmd_res['Atilde'].shape=}")
print(f"{dmd_res['u_r'].shape=}")
print(f"{dmd_res['L'].shape=}")
print(f"{dmd_res['W'].shape=}")
print(f"{dmd_res['Phi'].shape=}")
print(f"{dmd_res['amplitudes'][0].shape=}")

print('done')


A = dmd_res['A']
gramT = 3
gramPath = '/scratch/indikar_root/indikar0/jpic/subOptSS/2015/'

n = A.shape[0]
# Submodular optimization
S = []              # selected sensors
R = list(range(n))  # remaining sensors

At = {0: np.eye(A.shape[0])}
for t in range(1, gramT):
    At[t] = A @ At[t-1]

G = {}
for vx in R:
    print(str(vx) + '/' + str(n))
    C = getC(n, [vx]) # create C matrix
    C = sparse.csr_matrix(C)
    G[vx] = np.zeros_like(A) # construct new gramian
    for t in range(gramT):         # vary finite time
        G[vx] += (At[t].T @ C.T) @ (C @ At[t])
    df = pd.DataFrame(data=G[vx].astype(float))
    df.to_csv(gramPath + 'gramT_' + str(gramT) + '/' + str(vx) + '.csv', sep=',', header=False, float_format='%.4f', index=False)
    
    