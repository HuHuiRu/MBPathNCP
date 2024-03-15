import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score
import pandas as pd
from GIP.GIP_gate import GIP_kernel
from NCP.consistency_projection import NSP
import random
from sklearn.model_selection import StratifiedKFold
from WKNKN.wknkn import WKNKN

def prediction(CPI, ABP):   #用来得到预测矩阵
    a=0.95
    K=7
    eta=0.9
    Gip_cp_path=GIP_kernel(ABP)    #3427*3427
    Gip_path_cp=GIP_kernel(ABP.transpose())   #176*176
    ICP=(a*CPI+(1-a)*Gip_cp_path)

    new_ABP = WKNKN(ABP, ICP, Gip_path_cp, K, eta)

    # ICP=AB
    # ICP=Gip_cp_path
    # predicted_matrix = NSP(drug_similarity=2.4*ICP, atc_similarity=0.5*Gip_path_cp, adjacency_matrix=0.1*ABP).network_NSP()
    predicted_matrix = NSP(ICP_matrix=ICP, Gip_path_cp_matrix=Gip_path_cp,new_ABP_matrix=new_ABP).network_NSP()
    return predicted_matrix


def get_all_samples(conjunction):  # 输入：邻接矩阵
    pos = []  # 存放邻接矩阵为1的索引
    neg = []  # 存放邻接矩阵为0的索引
    for index in range(conjunction.shape[0]):
        for col in range(conjunction.shape[1]):
            if conjunction[index, col] == 1:
                pos.append([index, col, 1])
            else:
                neg.append([index, col, 0])
    pos_len = len(pos)

    # # 正负样本数量相同(负样本数量是正样本数量的1倍)
    # new_neg = random.sample(neg,pos_len)  # 随机获取和pos样本相同数量的neg样本
    # samples = pos + new_neg   #正负样本数量一样
    # neg_len = len(new_neg)

    # 所有的正负样本
    samples = pos + neg
    neg_len=len(neg)

    samples = random.sample(samples, len(samples))
    samples = np.array(samples)
    print("pos_len")
    print(pos_len)
    print("neg_len")
    print(neg_len)
    return samples



#预测CP_Pathway
def calculate_auc(CPI, CP_Pathway, predication):

    CPI = np.array(CPI)   #化合物和酶的关系矩阵
    CP_Pathway = np.array(CP_Pathway)  #化合物和酶与pathway的关联矩阵
    samples = get_all_samples(CP_Pathway)   #从CP_Pathway中获得同样数量的正负样本

    # 定义五倍交叉验证的折叠
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
    sum_auc_score = 0
    sum_average_precision=0
    all_auc_scores = []
    all_average_precision=[]

    # 进行交叉验证
    #train_index\test_index是训练集和验证集的索引
    #ABP_test_indices就是验证集
    labels=samples[:, 2]

    for train_index, test_index in kf.split(samples,labels):
        # train_samples, test_samples = samples[train_index], samples[test_index]
        train_labels, test_labels = labels[train_index], labels[test_index]


        test_indices = samples[test_index, :] #ABP_test_indices：(行号，列号，值)

        #令矩阵的训练部分为0
        CP_Pathway_train = CP_Pathway.copy()

        for i in test_indices:   #i 是什么？
            CP_Pathway_train[i[0], i[1]] = 0

        # 使用方法进行推荐矩阵预测   输入化合物和酶的关系矩阵，化合物和酶与Pathway的关系矩阵的训练部分（测试部分置0后的部分）
        predicted_matrix = predication(CPI, CP_Pathway_train)

        #获取测试集的索引
        test_sample_indices = test_indices[:, :2].astype(int)

        #从预测矩阵predicted_matrix中取出测试集部分（之前置为0的部分）

        predicted_matrix_test_labels = predicted_matrix[test_sample_indices[:, 0], test_sample_indices[:, 1]]

        #从原矩阵CP_Pathway中取出测试集部分

        CP_Pathway_test_labels = CP_Pathway[test_sample_indices[:, 0], test_sample_indices[:, 1]]


        #原矩阵的测试集与预测矩阵的测试集标签比较得出auc
        auc_score = roc_auc_score(CP_Pathway_test_labels, predicted_matrix_test_labels.flatten())
        average_precision = average_precision_score(CP_Pathway_test_labels, predicted_matrix_test_labels)


        # 输出AUC得分
        print(f"AUC score: {auc_score}")
        print("Average Precision:", average_precision)
        sum_auc_score += auc_score
        sum_average_precision+=average_precision
        all_auc_scores.append(auc_score)
        all_average_precision.append(average_precision)
        # 输出每个折叠的AUC得分
    for fold, score in enumerate(all_auc_scores):
        print(f"{score}")
    for fold, average_precision in enumerate(all_average_precision):
        print(f"{average_precision}")
    print(f"mean-AUC score: {sum_auc_score/10}")
    print(f"mean_average_precision score: {sum_average_precision / 10}")


#chemical_and_enzyme
CP_Pathway=pd.read_csv("Data/Pathways/Path_CPC_matrix.csv",index_col=0).to_numpy()
# CP_Pathway=pd.read_csv("Data/Pathways/CPC_Pathway.csv",index_col=0).to_numpy()
CPI=pd.read_csv('Data/chemical_enzyme_interactions/3427CPC_net.csv', index_col=0).to_numpy()

#only_chemical
# CP_Pathway=pd.read_csv("DataSet/CP_Pathway/Path_CID_matrix.csv",index_col=0).to_numpy() #cid_path
# CPI=pd.read_csv('DataSet/CPC_net/2329CCI_matrix.csv',index_col=0).to_numpy()  #cci



# #only_enzyme
# CP_Pathway=pd.read_csv("DataSet/CP_Pathway/Path_ENSP_matrix.csv",index_col=0).to_numpy()
# CPI=pd.read_csv('DataSet/CPC_net/1098PPI_matrix.csv',index_col=0).to_numpy()


# 计算AUC
calculate_auc(CPI, CP_Pathway, prediction)
