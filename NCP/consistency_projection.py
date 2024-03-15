import numpy as np


class NSP(object):
    def __init__(self, ICP_matrix, Gip_path_cp_matrix, new_ABP_matrix):
        self.ICP_matrix = ICP_matrix
        self.Gip_path_cp_matrix = Gip_path_cp_matrix
        self.new_ABP_matrix = new_ABP_matrix


    def PCP(self):
        temp_matrix = np.dot(self.new_ABP_matrix, self.Gip_path_cp_matrix)
        modulus = np.linalg.norm(self.new_ABP_matrix, axis=1).reshape(-1, 1)
        return temp_matrix / modulus

    def CECP(self):
        temp_matrix = np.dot(self.ICP_matrix, self.new_ABP_matrix)
        modulus = np.linalg.norm(self.new_ABP_matrix, axis=0).reshape(1, -1)
        return temp_matrix / modulus

    def calculate_modulus_sum(self):
        index_modulus = np.linalg.norm(self.ICP_matrix, axis=1).reshape(-1, 1)
        columns_modulus = np.linalg.norm(self.Gip_path_cp_matrix, axis=0).reshape(1, -1)
        return index_modulus + columns_modulus

    def network_NSP(self):
        result = np.nan_to_num((np.nan_to_num(self.PCP())+np.nan_to_num(self.CECP()))/self.calculate_modulus_sum())
        return result

