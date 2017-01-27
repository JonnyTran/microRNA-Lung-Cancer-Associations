import pandas as pd
import numpy as np
import networkx as nx
from networkx.algorithms import bipartite
from collections import OrderedDict
import operator

class miRNATargetNetwork:
    def __init__(self, dys_threshold=0.6):
        self.dys_threshold = dys_threshold
        self.B = nx.Graph()

    def add_miRNA_nodes(self, miRNAs):
        self.B.add_nodes_from(miRNAs, bipartite=0)

    def add_target_nodes(self, targets):
        self.B.add_nodes_from(targets, bipartite=1)

    def train(self, miRNAs_tumor, targets_tumor, miRNAs_normal, targets_normal, putative_assocs):
        """
        Constructing the MTDN from xu2011prioritizing

        :param miRNAs_tumor: array of shape (n_samples, n_miRNAs) for tumor samples
        :param targets_tumor: array of shape (n_samples, n_genes) for tumor samples
        :param miRNAs_normal: array of shape (n_samples, n_miRNAs) for normal  samples
        :param targets_normal: array of shape (n_samples, n_genes) for normal samples
        """
        miRNAs = miRNAs_tumor.columns
        targets = targets_normal.columns
        self.add_miRNA_nodes(miRNAs)
        self.add_target_nodes(targets)

        n_A = miRNAs_tumor.shape[0]
        n_B = miRNAs_normal.shape[0]
        print 'n_A', n_A
        print 'n_B', n_B

        for i in putative_assocs.index:
            m = putative_assocs.ix[i]['MiRBase ID']
            t = putative_assocs.ix[i]['Gene Symbol']

            if (m in miRNAs) and (t in targets):
                miRNA_target_A_corr = np.dot(miRNAs_tumor[m] - np.mean(miRNAs_tumor[m]),
                                             targets_tumor[t] - np.mean(targets_tumor[t])) / \
                                      ((n_A - 1) * np.std(miRNAs_tumor[m]) * np.std(targets_tumor[t]))

                miRNA_target_B_corr = np.dot(miRNAs_normal[m] - np.mean(miRNAs_normal[m]),
                                             targets_normal[t] - np.mean(targets_normal[t])) / \
                                      ((n_B - 1) * np.std(miRNAs_normal[m]) * np.std(targets_normal[t]))

                dys = miRNA_target_A_corr - miRNA_target_B_corr
                # print m, '<->', t, ':', dys

                if abs(dys) >= self.dys_threshold:
                    self.B.add_edge(m, t, dys=dys)

                    # if miRNA_target_A_corr <= -self.dys_threshold:
                    #     self.B.add_edge(m, t, tumor_corr=miRNA_target_A_corr)
                    #
                    # if miRNA_target_B_corr <= -self.dys_threshold:
                    #     self.B.add_edge(m, t, tumor_corr=miRNA_target_B_corr)

                    # for m in miRNAs:
                    #     edge_count = 1
                    #     miRNA_A_m = miRNAs_A[m] - np.mean(miRNAs_A[m])
                    #     miRNA_B_m = miRNAs_B[m] - np.mean(miRNAs_B[m])
                    #     miRNA_A_m_std = np.std(miRNAs_A[m])
                    #     miRNA_B_m_std = np.std(miRNAs_B[m])
                    #     for t in targets:
                    #         miRNA_target_A_corr = np.dot(miRNA_A_m, targets_A[t] - np.mean(targets_A[t])) / \
                    #                               ((n_A - 1) * miRNA_A_m_std * np.std(targets_A[t]))
                    #
                    #         miRNA_target_B_corr = np.dot(miRNA_B_m, targets_B[t] - np.mean(targets_B[t])) / \
                    #                               ((n_B - 1) * miRNA_B_m_std * np.std(targets_B[t]))
                    #
                    #         dys = miRNA_target_A_corr - miRNA_target_B_corr
                    #         # print m, '<->', t, ':', dys
                    #
                    #         if abs(dys) >= self.threshold:
                    #             self.add_edge(m, t, dys=dys)
                    #             edge_count += 1
                    #
                    #     print m, ':', edge_count

    def get_miRNA_groups(self, mirna_list, smaller_groups=True):
        miRNAs_nodes = set(n for n, d in self.B.nodes(data=True) if d['bipartite'] == 0)
        targets_nodes = set(self.B) - miRNAs_nodes

        targets_nodes_degrees = nx.bipartite.degrees(self.B, targets_nodes)[1]

        sorted_targets_nodes_degrees = sorted(targets_nodes_degrees.items(), key=operator.itemgetter(1),
                                              reverse=smaller_groups)

        mirna_group_assg = OrderedDict((miRNA, -1) for miRNA in mirna_list)
        self.miRNA_groups = []
        self.miRNA_groups_int = []

        # For every target, find its neighbor miRNA's, which forms a group.
        # The miRNA's in a groups are then assigned a corresponding number
        group_counter = 1
        for (target, n_neighbors) in sorted_targets_nodes_degrees:
            if n_neighbors > 1:
                target_neighbors = self.B.neighbors(target)
                self.miRNA_groups.append(target_neighbors)
                self.miRNA_groups_int.append([mirna_list.index(miRNA) for miRNA in target_neighbors])

                for miRNA in target_neighbors:
                    mirna_group_assg[miRNA] = group_counter
                group_counter += 1

        groups_unique = np.unique(mirna_group_assg.values())
        group_counter = 1

        # Ensure the rest of the miRNAs not in groups to have unique group number, starting from 1
        for miRNA, group_assg in mirna_group_assg.iteritems():
            if group_assg == -1:
                while group_counter in groups_unique:
                    group_counter += 1
                mirna_group_assg[miRNA] = group_counter
                group_counter += 1

        return mirna_group_assg.values()

    def find_miRNA_groups(self, miRNA):
        groups = []
        for group in self.miRNA_groups:
            if miRNA in group:
                groups.append(len(group))

        return groups
