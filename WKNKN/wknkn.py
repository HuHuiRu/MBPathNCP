import numpy as np
import pandas as pd
from GIP.GIP_gate import GIP_kernel
def WKNKN(A, SD, ST, K, eta):
    Ad = np.zeros(A.shape)
    At = np.zeros(A.shape)
    wi = np.zeros((K,))
    wj = np.zeros((K,))
    num_ce, num_p = A.shape
    for i in np.arange(num_ce):
        dnn_i = np.argsort(SD[i,:])[::-1][1:K+1]
        Zd = np.sum(SD[i, dnn_i])
        for ii in np.arange(K):
            wi[ii] = (eta ** (ii)) * SD[i,dnn_i[ii]]
        if not np.isclose(Zd, 0.):
            Ad[i,:] = np.sum(np.multiply(wi.reshape((K,1)), A[dnn_i,:]), axis=0) / Zd
    for j in np.arange(num_p):
        tnn_j = np.argsort(ST[j, :])[::-1][1:K+1]
        Zt = np.sum(ST[j, tnn_j])
        for jj in np.arange(K):
            wj[jj] = (eta ** (jj)) * ST[j,tnn_j[jj]]
        if not np.isclose(Zt, 0.):
            At[:,j] = np.sum(np.multiply(wj.reshape((1,K)), A[:,tnn_j]), axis=1) / Zt
    Adt = (Ad + At)/2
    x, y = np.where(Adt > A)

    A_tem = A.copy()
    A_tem[x, y] = Adt[x, y]
    return A_tem


