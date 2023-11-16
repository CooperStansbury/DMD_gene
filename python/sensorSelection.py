import pandas as pd
import numpy as np
from copy import deepcopy
import os
import sys
import scipy

import numpy as np

def submodularSensorSelection(A, gramT=1,maxSensors=2, subCriteria):
    n = A.shape[0]
    # Submodular optimization
    S = []              # selected sensors
    R = list(range(n))  # remaining sensors

    # while selecting more sensors
    while len(S) < maxSensors:
        M = np.zeros(len(R))  # save scores for each sensor
        # try each of the remaining sensors
        for i, vx in enumerate(R):
            C = getC(n, np.append(S, vx))  # create C matrix
            G = np.zeros_like(A)           # construct new gramian
            for t in range(gramT):         # vary finite time
                G += np.dot(np.dot(A.T, C.T), np.dot(C, A))
            if subCriteria == 1:
                M[i] = np.trace(G)      # Four measures of submodularity
            elif subCriteria == 2:
                M[i] = np.trace(np.linalg.inv(G))
            elif subCriteria == 3:
                M[i] = np.log(np.linalg.det(G))
            elif subCriteria == 4:
                M[i] = np.linalg.matrix_rank(G)
        vx = np.argmax(M)  # find highest weighted next sensor
        S.append(R[vx])    # select the next sensor
        R.pop(vx)          # remove sensor from remaining vertices
    return S

def getC(n, idxs):
    # Define your getC function if it's not already defined
    C = sp.sparse(n, len(idxs))
    for i in range(idxs):
        C[i, idxs[i]] = 1
    return C


