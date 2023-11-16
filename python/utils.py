import os
import sys
import pandas as pd
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import scipy
# import scanpy as sc
from scipy import stats
from Bio import SeqIO
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
import leidenalg
import io


def cluster(X, knn=11, resolution=1.0):
    """A function to cluster a data matrix 

    Params:
    --------------
    X (np.array):
        The data matrix
    knn (int):
        The number of nearest neighbors
    resolution (float):
        The resolution parameter for the Leiden clustering 
        algorithm
        
    Returns:
    --------------
    groups (list):
        Cluster membership Ids
    
    """
    # transform to adjacency graph 
    nbrs = NearestNeighbors(n_neighbors=knn, algorithm='ball_tree').fit(X)
    g = nbrs.kneighbors_graph(X).toarray()
    g = sc._utils.get_igraph_from_adjacency(g, directed=False)
    
    partition_type = leidenalg.RBConfigurationVertexPartition
    part = leidenalg.find_partition(g, partition_type)
    groups = np.array(part.membership)

    return groups


def ncolor(n, cmap='viridis'):
    """A function to generate a list of hex codes from a colormap 

    Params:
    --------------
    n (int):
        The number of colors
    cmap (str):
        A colormap string
        
    Returns:
    --------------
    colors (list):
        A list of hex codes
    """
    cmap = matplotlib.cm.get_cmap(cmap)
    arr = np.linspace(0, 1, n)
    return [matplotlib.colors.rgb2hex(cmap(x)) for x in arr] 


def parseKEGG(pathId):
    """A function to return a list of gene names given a KEGG
    pathway id
    
    Params:
    --------------
    pathId (str):
        A pathway id, e.g., 'hsa04260' for Cardiac muscle contraction - Homo sapiens (human)

    Returns:
    --------------
    genes (list):
        A list of genes from the pathway
    """
    genes = []
    results = REST.kegg_get(pathId).read()
    current_section = None
    for line in results.rstrip().split("\n"):
        section = line[:12].strip()  # section names are within 12 columns
        if not section == "":
            current_section = section

        if current_section == "GENE":
            linesplit = line[12:].split("; ")
            gene_identifiers = linesplit[0]
            gene_id, gene_symbol = gene_identifiers.split()
    
            if not gene_symbol in genes:
                genes.append(gene_symbol)
    return genes


def getPathname(pathId):
    """A function to return the KEGG pathway name
    
    Params:
    --------------
    pathId (str):
        A pathway id, e.g., 'hsa04260' for Cardiac muscle contraction - Homo sapiens (human)

    Returns:
    --------------
    pathname (str):
        The name of the pathway as a string
    """
    result = REST.kegg_list(pathId).read()
    return result.split("\t")[1].split("-")[0].strip()


def makeColorbar(cmap, width, hieght, title, orientation, tickLabels):
    """A function to build a colorbar.

    Params:
    --------------
    cmap (str):
        A colomap map string label 
    width (float):
        The width of the colorbar
    hieght (float):
        The hieght of the colorbar
    title (str):
        The title of the colorbar
    orientation (str):
        Either 'vertical' or 'horizontal'
    tickLabels (list of str):
        Evenly spaced tick labels
    
    Returns:
    --------------
    colobar (plt.Figure):
        A colorbar image
    """
    a = np.array([[0,1]])
    plt.figure(figsize=(width, hieght))
    img = plt.imshow(a, cmap=cmap)
    plt.gca().set_visible(False)
    cax = plt.axes([0.1, 0.2, 0.8, 0.6])
    ticks = np.linspace(0,1 , len(tickLabels))
    cbar = plt.colorbar(orientation=orientation, 
                        cax=cax, 
                        label=title,
                        ticks=ticks)

    if orientation == 'vertical':
        cbar.ax.set_yticklabels(tickLabels)
    else:
        cbar.ax.set_xticklabels(tickLabels)

        
def _normalize_data(X, counts, after=None, copy=False):
    """A function to normalize data: adapted from scanpy """
    X = X.copy() if copy else X
    if issubclass(X.dtype.type, (int, np.integer)):
        X = X.astype(np.float32)  # TODO: Check if float64 should be used
    else:
        counts_greater_than_zero = counts[counts > 0]

    after = np.median(counts_greater_than_zero, axis=0) if after is None else after
    counts += counts == 0
    counts = counts / after
    if scipy.sparse.issparse(X):
        sparsefuncs.inplace_row_scale(X, 1 / counts)
    elif isinstance(counts, np.ndarray):
        np.divide(X, counts[:, None], out=X)
    else:
        X = np.divide(X, counts[:, None])  # dask does not support kwarg "out"
    return X


def normalize(df, target_sum=1):
    """A function to normalize data: adapted from scanpy  """
    index = df.index
    columns = df.columns
    X = df.to_numpy().copy()
    counts_per_cell = X.sum(1)
    counts_per_cell = np.ravel(counts_per_cell)
    cell_subset = counts_per_cell > 0
    Xnorm = _normalize_data(X, counts_per_cell, target_sum)
    
    ndf = pd.DataFrame(Xnorm, columns=columns, index=index)
    return ndf