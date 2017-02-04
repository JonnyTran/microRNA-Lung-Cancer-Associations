import pandas
import numpy as np
import networkx as nx
from networkx.algorithms import bipartite
from collections import OrderedDict
import operator
import pandas
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool
from sklearn.cluster import AgglomerativeClustering
from multiprocessing import Manager

class miRNATargetNetwork:
    def __init__(self, miRNAs, targets, dys_threshold=0.6):
        self.dys_threshold = dys_threshold
        self.B = nx.DiGraph()
        self.mirna_list = miRNAs
        self.genes_list = targets
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
        n_A = miRNA_A.shape[0]
        n_B = miRNA_B.shape[0]
        print 'n_A', n_A
        print 'n_B', n_B

        # mgr = Manager()
        # ns = mgr.Namespace()
        # ns.putative_assocs = putative_assocs
        # ns.dys_threshold = self.dys_threshold
        # ns.miRNA_A = miRNA_A
        # ns.gene_A = gene_A
        # ns.miRNA_B = miRNA_B
        # ns.gene_B = gene_B
        # ns.n_A = n_A
        # ns.n_B = n_B
        # ns.tag = tag
        #
        # def test(i):
        #     m = ns.putative_assocs.ix[i]['MiRBase ID']
        #     t = ns.putative_assocs.ix[i]['Gene Symbol']
        #
        #     miRNA_gene_A_corr = np.dot(ns.miRNA_A[m] - np.mean(ns.miRNA_A[m]),
        #                                ns.gene_A[t] - np.mean(ns.gene_A[t])) / \
        #                         ((n_A - 1) * np.std(ns.miRNA_A[m]) * np.std(ns.gene_A[t]))
        #
        #     miRNA_gene_B_corr = np.dot(ns.miRNA_B[m] - np.mean(ns.miRNA_B[m]),
        #                                ns.gene_B[t] - np.mean(ns.gene_B[t])) / \
        #                         ((n_B - 1) * np.std(ns.miRNA_B[m]) * np.std(ns.gene_B[t]))
        #
        #     dys = miRNA_gene_A_corr - miRNA_gene_B_corr
        #     # print m, '<->', t, ':', dys
        #
        #     if abs(dys) >= ns.dys_threshold:
        #         return (m, t, dys, ns.tag)
        #
        # pool = Pool(8)
        # results = pool.map(test, putative_assocs.index.tolist())
        #
        # edges_added = 0
        # for tup in results:
        #     if tup:
        #         edges_added += 1
        #         print tup
        #         self.B.add_edge(tup[0], tup[1], dys=tup[2], tag=tup[3])


        edges_added = 0
        for row in putative_assocs.iterrows():
            m = row[1]['MiRBase ID']
            t = row[1]['Gene Symbol']

            miRNA_gene_A_corr = np.dot(miRNA_A[m] - np.mean(miRNA_A[m]),
                                       gene_A[t] - np.mean(gene_A[t])) / \
                                ((n_A - 1) * np.std(miRNA_A[m]) * np.std(gene_A[t]))

            miRNA_gene_B_corr = np.dot(miRNA_B[m] - np.mean(miRNA_B[m]),
                                       gene_B[t] - np.mean(gene_B[t])) / \
                                ((n_B - 1) * np.std(miRNA_B[m]) * np.std(gene_B[t]))

            dys = miRNA_gene_A_corr - miRNA_gene_B_corr
            # print m, '<->', t, ':', dys

            if abs(miRNA_gene_A_corr) > 1.0 or abs(miRNA_gene_B_corr) > 1.0:
                print 'we got a problem', miRNA_gene_A_corr, miRNA_gene_B_corr

            if abs(dys) >= self.dys_threshold:
                self.B.add_edge(m, t, dys=dys, tag=tag)
                edges_added += 1

        return edges_added

    def build_miRNA_features(self):
        tags = ['normal-StgI', 'StgI-StgII', 'StgII-StgIII', 'StgIII-StgIV']

        normal_stgI_genes = np.unique(
            [tup[1] for tup in self.B.edges(data=True) if tup[2]['tag'] == 'normal-StgI']).tolist()
        stgI_StgII_genes = np.unique(
            [tup[1] for tup in self.B.edges(data=True) if tup[2]['tag'] == 'StgI-StgII']).tolist()
        stgII_StgIII_genes = np.unique(
            [tup[1] for tup in self.B.edges(data=True) if tup[2]['tag'] == 'StgII-StgIII']).tolist()
        stgIII_StgIV_genes = np.unique(
            [tup[1] for tup in self.B.edges(data=True) if tup[2]['tag'] == 'StgIII-StgIV']).tolist()
        miRNAs_in_MTDN = np.unique([tup[0] for tup in self.B.edges()])

        print 'normal_stgI_genes', len(normal_stgI_genes)
        print 'stgI_StgII_genes', len(stgI_StgII_genes)
        print 'stgII_StgIII_genes', len(stgII_StgIII_genes)
        print 'stgIII_StgIV_genes', len(stgIII_StgIV_genes)
        print 'miRNAs_in_MTDN', len(miRNAs_in_MTDN)

        normal_stgI_genes_names = [gene + '/normal-StgI' for gene in normal_stgI_genes]
        stgI_StgII_genes_names = [gene + '/StgI-StgII' for gene in stgI_StgII_genes]
        stgII_StgIII_genes = [gene + '/StgII-StgIII' for gene in stgII_StgIII_genes]
        stgIII_StgIV_genes = [gene + '/StgIII-StgIV' for gene in stgIII_StgIV_genes]

        self.miRNA_target_assn_matrix = pandas.DataFrame(columns=normal_stgI_genes_names +
                                                                 stgI_StgII_genes_names +
                                                                 stgII_StgIII_genes +
                                                                 stgIII_StgIV_genes)

        for miRNA in miRNAs_in_MTDN:
            self.miRNA_target_assn_matrix.loc[miRNA] = 0

        for edge in self.B.edges(data=True):
            self.miRNA_target_assn_matrix.loc[edge[0]][edge[1] + '/' + edge[2]['tag']] = 1

    def run_miRNA_clustering(self, n_cluster=20, linkage='complete'):

        mirna_cluster = AgglomerativeClustering(n_clusters=n_cluster, affinity='l1', linkage=linkage).fit(
            self.miRNA_target_assn_matrix)
        self.miRNA_cluster_assgn = mirna_cluster.fit_predict(self.miRNA_target_assn_matrix)

        self.miRNA_clusters_int = []
        for cluster_idx in range(n_cluster):
            self.miRNA_clusters_int.append([])
            for mirna in [self.miRNA_target_assn_matrix.index[mirna_idx]
                          for mirna_idx, cluster_assg in enumerate(self.miRNA_cluster_assgn) if
                          cluster_assg == cluster_idx]:
                self.miRNA_clusters_int[cluster_idx].append(self.mirna_list.index(mirna))

        return np.bincount(mirna_cluster.fit_predict(self.miRNA_target_assn_matrix))

    def get_miRNA_cluster_assgn(self):
        mirna_group_assg = OrderedDict((miRNA, -1) for miRNA in self.mirna_list)

        for mirna, cluster_assgn in zip(self.miRNA_target_assn_matrix.index,
                                        self.miRNA_cluster_assgn):
            mirna_group_assg[mirna] = cluster_assgn

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
