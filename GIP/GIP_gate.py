import numpy as np
import pandas as pd
def GIP_kernel(Asso_RNA_Dis):
    nc = Asso_RNA_Dis.shape[0]
    # initate a matrix as results matrix 初始化一个矩阵作为结果矩阵
    matrix = np.zeros((nc, nc))
    # calculate the down part of GIP fmulate   #计算GIP公式的下部分
    r = getGosiR(Asso_RNA_Dis)
    # calculate the results matrix   计算结果矩阵
    for i in range(nc):
        for j in range(nc):
            # calculate the up part of GIP formulate  计算GIP公式的上部分
            temp_up = np.square(np.linalg.norm(Asso_RNA_Dis[i, :] - Asso_RNA_Dis[j, :]))
            if r == 0:
                matrix[i][j] = 0
            elif i == j:
                matrix[i][j] = 1
            else:
                matrix[i][j] = np.e ** (-temp_up / r)
    return matrix

def getGosiR(Asso_RNA_Dis):
    # calculate the r in GOsi Kerel
    nc = Asso_RNA_Dis.shape[0]
    summ = 0
    for i in range(nc):
        x_norm = np.linalg.norm(Asso_RNA_Dis[i, :])
        x_norm = np.square(x_norm)
        summ = summ + x_norm
    r = summ / nc
    return r
