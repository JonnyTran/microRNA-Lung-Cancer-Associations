import pandas as pd
import numpy as np
import networkx as nx
from networkx.algorithms import bipartite


class miRNATargetNetwork:
    def __init__(self, threshold=0.6):
        self.threshold = threshold
        self.B = nx.Graph()

    def add_miRNA_nodes(self, miRNAs):
        self.B.add_nodes_from(miRNAs, bipartite=0)

    def add_target_nodes(self, targets):
        self.B.add_nodes_from(targets, bipartite=1)

    def train(self, miRNAs_A, targets_A, miRNAs_B, targets_B):
        # Constructing the MTDN from xu2011prioritizing
        """

        :param miRNAs_A: Pandas dataframe
        :param targets_A: Pandas dataframe
        :param miRNAs_B: Pandas dataframe
        :param targets_B: Pandas dataframe
        """
        miRNAs = miRNAs_A.columns
        targets = targets_A.columns
        self.add_miRNA_nodes(miRNAs)
        self.add_target_nodes(targets)

        n_A = miRNAs_A.shape[0]
        n_B = miRNAs_B.shape[0]

        for m in miRNAs:
            for t in targets:
                dys = np.dot((miRNAs_A[m] - np.mean(miRNAs_A[m])), (targets_A[t] - np.mean(targets_A[t]))) / (
                (n_A - 1) * np.std(miRNAs_A[m]) * np.std(targets_A[t])) - \
                      np.dot((miRNAs_B[m] - np.mean(miRNAs_B[m])), (targets_B[t] - np.mean(targets_B[t]))) / (
                      (n_B - 1) * np.std(miRNAs_B[m]) * np.std(targets_B[t]))
                print m, '-', t, ':', dys
                if dys >= self.threshold:
                    self.add_edge(m, t, dys=dys)

    def add_edge(self, miRNA, target, dys):
        self.B.add_edge(miRNA, target, dys=dys)
