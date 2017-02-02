import pandas as pd
import numpy as np
import networkx as nx
from networkx.algorithms import bipartite
from collections import OrderedDict
import operator
from multiprocessing import Pool


class miRNATargetNetwork:
    def __init__(self, miRNAs, targets, dys_threshold=0.6):
        self.dys_threshold = dys_threshold
        self.B = nx.Graph()
        self.mirna_list = miRNAs.tolist()
        self.add_miRNA_nodes(miRNAs)
        self.add_target_nodes(targets)

    def add_miRNA_nodes(self, miRNAs):
        self.B.add_nodes_from(miRNAs, bipartite=0)

    def add_target_nodes(self, targets):
        self.B.add_nodes_from(targets, bipartite=1)

    def fit(self, miRNA_A, gene_A, miRNA_B, gene_B, putative_assocs, tag=""):
        """
        Constructing the MTDN from xu2011prioritizing

        :param miRNA_A: array of shape (n_samples, n_miRNAs) for tumor samples
        :param gene_A: array of shape (n_samples, n_genes) for tumor samples
        :param miRNA_B: array of shape (n_samples, n_miRNAs) for normal  samples
        :param gene_B: array of shape (n_samples, n_genes) for normal samples
        """
        miRNAs = miRNA_A.columns
        targets = gene_B.columns

        n_A = miRNA_A.shape[0]
        n_B = miRNA_B.shape[0]
        print 'n_A', n_A
        print 'n_B', n_B

        edges_added = 0

        for i in putative_assocs.index:
            m = putative_assocs.ix[i]['MiRBase ID']
            t = putative_assocs.ix[i]['Gene Symbol']

            if (m in miRNAs) and (t in targets):
                miRNA_gene_A_corr = np.dot(miRNA_A[m] - np.mean(miRNA_A[m]),
                                           gene_A[t] - np.mean(gene_A[t])) / \
                                    ((n_A - 1) * np.std(miRNA_A[m]) * np.std(gene_A[t]))

                miRNA_gene_B_corr = np.dot(miRNA_B[m] - np.mean(miRNA_B[m]),
                                           gene_B[t] - np.mean(gene_B[t])) / \
                                    ((n_B - 1) * np.std(miRNA_B[m]) * np.std(gene_B[t]))

                dys = miRNA_gene_A_corr - miRNA_gene_B_corr
                # print m, '<->', t, ':', dys

                if abs(dys) >= self.dys_threshold:
                    self.B.add_edge(m, t, dys=dys, tag=tag)
                    edges_added += 1

        return edges_added

    def get_miRNA_features(self):
        miRNAs_nodes = set(n for n, d in self.B.nodes(data=True) if d['bipartite'] == 0)
        targets_nodes = set(self.B) - miRNAs_nodes

        targets_nodes_degrees = nx.bipartite.degrees(self.B, targets_nodes)[1]

        edges = self.B.edges()

    def get_miRNA_group_assgn(self, smaller_groups=True):
        miRNAs_nodes = set(n for n, d in self.B.nodes(data=True) if d['bipartite'] == 0)
        targets_nodes = set(self.B) - miRNAs_nodes

        targets_nodes_degrees = nx.bipartite.degrees(self.B, targets_nodes)[1]

        sorted_targets_nodes_degrees = sorted(targets_nodes_degrees.items(), key=operator.itemgetter(1),
                                              reverse=smaller_groups)

        mirna_group_assg = OrderedDict((miRNA, -1) for miRNA in self.mirna_list)
        self.miRNA_groups = []
        self.miRNA_groups_int = []

        # For every target, find its neighbor miRNA's, which forms a group.
        # The miRNA's in a groups are then assigned a corresponding number
        group_counter = 1
        for (target, n_neighbors) in sorted_targets_nodes_degrees:
            if n_neighbors > 1:
                target_neighbors = self.B.neighbors(target)
                self.miRNA_groups.append(target_neighbors)
                self.miRNA_groups_int.append([self.mirna_list.index(miRNA) for miRNA in target_neighbors])

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
        for i, group in enumerate(self.miRNA_groups):
            if miRNA in group:
                groups.append(i)

        return groups
