import pandas as pd
import numpy as np
from copy import deepcopy
import os
import sys


def data2dmd(data):
    """A function to stack replicates for DMD as in Hasnain et al.

    Params:
    --------------
    data (pd.DataFrame):
        A dataframe with genes as columns and rows as timepoints * replicates

    Returns:
    --------------
    dmd_data (np.array):
        An array of shape (genes, timepoints, replicates)
    """
    labels = pd.DataFrame(data.index, columns=['expId'])
    labels['time'] = labels['expId'].apply(lambda x: x.split("_")[0])
    labels['replicate'] = labels['expId'].apply(lambda x: x.split("_")[1])

    replicates = sorted(labels['replicate'].unique())

    dmd_data = []

    for r in replicates:
        # get the indices of the replicate
        ind = labels[labels['replicate'] == r]['expId'].to_list()

        t_data = data.loc[ind]
        dmd_data.append(t_data)
        
    dmd_data = np.asarray(dmd_data)
    dmd_data = np.swapaxes(dmd_data, 0, 2)
    return dmd_data
    

def getOHT(u, s, vh):
    """A function to compute the optimal hard threshold from the SVD 
    of the data. NOTE: assumes a tall, skinny matrix 

    Params:
    --------------
    u (np.array):
        Left singular vectors
    s (np.array):
        Diagonal matrix of singular values
    vh (np.array):
        Right singular vectors transposed

    Returns:
    --------------
    oht (int):
        The index of the optimal hardthreshold from a non-square
        martrix with noise level unknown.
    """
    n = u.shape[0]
    m = vh.shape[0] 
    
    beta = m / n
    omega = (0.56*beta**3) - (0.95 * beta**2) + (1.82 * beta) + 1.43
    y_med = np.median(s)
    tau = omega * y_med
    s_ind = np.argwhere(s >= tau)
    oht = np.max(s_ind) 
    return oht
    

def dmd_reshape(data):
    """A utility function to reshape the data as in hasnain el at.
     
    Params:
    --------------
    data (np.array):
        An array of shape (genes, timepoints, replicates)

    Returns:
    --------------
    Xp (np.array):
        The first m-1 timepoints for all replicates

    Xf (np.array):
        The last m-1 timepoints for all replicates
    """
    n, m, r = data.shape

    Xp = data[:,:-1].reshape(n, (m-1)*r, order='F') 
    Xf = data[:,1:].reshape(n, (m-1)*r, order='F') 

    return Xp, Xf


def n_step_prediction(A, X, ntimepts, nreps):
    """A function from  on Hasnain et al. 2023 to predict using the DMD solution
    
    Params:
    --------------
    A (np.array):
        The learned DMD operator
    X (np.array): 
        the data 
    ntimepts (int):
        The number of prediction steps
    nreps (int):
        The number of replicates

    
    Returns:
    --------------  
    X_pred (np.array):
        The DMD predicted matrix
    cd (float):
        The correlation of the data with the prediction
    """
    X = X[:,:ntimepts].reshape(len(X), (ntimepts)*nreps,order='F')
    X_pred = np.zeros((A.shape[0],ntimepts*nreps))
    count = 0
    for i in range(0,nreps):
        x_test_ic = X[:,i*(ntimepts):i*(ntimepts)+1]
        for j in range(0,ntimepts):
            X_pred[:,count:count+1] = np.dot(np.linalg.matrix_power(A,j),x_test_ic) 
            count += 1
    feature_means = np.mean(X,axis=1).reshape(len(X),1)
    cd = 1 - ((np.linalg.norm(X - X_pred,ord=2)**2)/(np.linalg.norm(X - feature_means,ord=2)**2))   # coeff of determination aka R^2 
    return X_pred, cd


def embed_data(data, u, rank):
    """A utility function to embed the data based on the 
    low-rank approximation of Xp 
    
    Params:
    --------------
    data (np.array):
        An array of shape (genes, timepoints, replicates)
    u (np.array):
        The left singular vectors of Xp
    rank (int):
        The rank truncation for u
    
    Returns:
    --------------  
    data_embedded (np.array):
        The embedded data
    """
    u_r = u[:, 0:rank] # truncate to rank-r
    n, m, r = data.shape
    
    data_embedded = np.zeros((rank, m, r))

    for i in range(r):
        data_embedded[:,:,i] = np.dot(u_r.T, data[:,:,i])
    return data_embedded



def dmd(data, exact=False, rank=None):
    """A function to compute the DMD of the data based on Hasnain et al. 2023
    
    Params:
    --------------
    data (np.array):
        An array of shape (genes, timepoints, replicates)
    exact (bool):
        If `True' rank is the minimum of (ncol, nrow) and the DMD solution is
        the regression of the data. If `False', expects a rank paramater
    rank (int or None):
        If `None', the trucated SVD will be computed using the optimal hard threshold
    
    Returns:
    --------------
    dmd_res (dict):
        A dictionary of DMD results
    """
    n, m, r = data.shape
    # RESHAPE DATA
    Xp, Xf = dmd_reshape(data)

    # perform DMD
    if exact:
        A = Xf @ np.linalg.pinv(Xp)  # full model A
        rank = np.minimum(Xp.shape[1], Xp.shape[0])
    else:
        u, s, vh = np.linalg.svd(Xp)

        if rank == None: 
            rank = getOHT(u, s, vh)

        # perform DMD
        u_r = u[:, 0:rank] # truncate to rank-r
        s_r = s[0:rank]
        vh_r = vh[0:rank, :]
        Atilde = u_r.T @ Xf @ vh_r.T @ np.diag(1/s_r) # low-rank dynamics
        A = u_r@Atilde@u_r.T
        
    # DMD EIGENVALUES AND EIGENVECTORS
    L, W = np.linalg.eig(Atilde)
    
    # DMD MODES
    Phi = Xf @ vh_r.T @ np.diag(1/s_r) @ W
    Phi_hat = np.dot(u_r, W)

    # DMD AMPLITUDES: X0 in the eigenvector (of A) basis
    data_embedded = embed_data(data, u, rank)
    amps = []
    for i in range(r):
        b_ri = np.linalg.inv(np.dot(W, np.diag(L))) @ data_embedded[:,:,i]
        amps.append(b_ri)
    
    return {
        'A' : A,
        'Atilde' : Atilde,
        'rank' : rank,
        'u_r' : u_r,
        'SVD' : [u, s, vh],
        'L' : L,
        'W' : W,
        'data_embedded' : data_embedded,
        'Phi' : Phi,
        'Phi_hat' : Phi_hat,
        'amplitudes' : amps
    }