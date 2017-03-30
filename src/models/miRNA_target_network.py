import copy
import operator
from collections import OrderedDict

import community  # python-louvain
import dask.dataframe as dd
import networkx as nx
import numpy as np
import pandas
import scipy.stats
from dask.multiprocessing import get
from sklearn.cluster import AgglomerativeClustering


class miRNATargetNetwork:
    def __init__(self, miRNAs, targets):
        self.B = nx.DiGraph()
        self.mirna_list = miRNAs
        self.genes_list = targets
        self.add_miRNA_nodes(miRNAs)
        self.add_target_nodes(targets)

    def add_miRNA_nodes(self, miRNAs):
        self.B.add_nodes_from(miRNAs, bipartite=0)

    def add_target_nodes(self, targets):
        self.B.add_nodes_from(targets, bipartite=1)

    def fit(self, miRNA_A, gene_A, miRNA_B, gene_B, putative_assocs, p_threshold=0.001, tag="A-B", n_jobs=5):
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

        edges_added = 0
        if putative_assocs is not None:
            putative_dd = dd.from_pandas(putative_assocs, npartitions=n_jobs)

            def calc_dys_A_B(df):
                result = []
                for row in df.iterrows():
                    m = row[1]['MiRBase ID']
                    t = row[1]['Gene Symbol']
                    miRNA_gene_A_corr = np.dot(miRNA_A[m] - np.mean(miRNA_A[m]),
                                               gene_A[t] - np.mean(gene_A[t])) / \
                                        ((n_A - 1) * np.std(miRNA_A[m]) * np.std(gene_A[t]))
                    miRNA_gene_B_corr = np.dot(miRNA_B[m] - np.mean(miRNA_B[m]),
                                               gene_B[t] - np.mean(gene_B[t])) / \
                                        ((n_B - 1) * np.std(miRNA_B[m]) * np.std(gene_B[t]))
                    dys = miRNA_gene_A_corr - miRNA_gene_B_corr
                    p_value = self.z_to_p_value(self.fisher_r_to_z(miRNA_gene_A_corr, n_A, miRNA_gene_B_corr, n_B))

                    if p_value <= p_threshold:
                        result.append((m, t, p_value))
                return result

            res = putative_dd.map_partitions(calc_dys_A_B, meta=putative_dd).compute(get=get)

            for res_partition in res:
                for tup in res_partition:
                    self.B.add_edge(tup[0], tup[1], dys=tup[2], tag=tag)
                    edges_added += 1
        ## Iterate through every miRNA-gene associations
        else:
            for m in self.mirna_list:
                m_A = miRNA_A[m] - np.mean(miRNA_A[m])
                m_B = miRNA_B[m] - np.mean(miRNA_B[m])
                std_A = np.std(miRNA_A[m])
                std_B = np.std(miRNA_B[m])
                for t in self.genes_list:
                    miRNA_gene_A_corr = np.dot(m_A,
                                               gene_A[t] - np.mean(gene_A[t])) / \
                                        ((n_A - 1) * std_A * np.std(gene_A[t]))

                    miRNA_gene_B_corr = np.dot(m_B,
                                               gene_B[t] - np.mean(gene_B[t])) / \
                                        ((n_B - 1) * std_B * np.std(gene_B[t]))

                    dys = miRNA_gene_A_corr - miRNA_gene_B_corr

                    if abs(dys) >= p_threshold and (abs(miRNA_gene_A_corr) > 0.3 or abs(miRNA_gene_B_corr) > 0.3):
                        self.B.add_edge(m, t, dys=dys, tag=tag)
                        edges_added += 1

        return edges_added

    def fisher_r_to_z(self, r_A, n_A, r_B, n_B):
        r_A_plus = 1.0 * r_A + 1
        r_A_minus = 1.0 - r_A
        r_B_plus = 1.0 * r_B + 1
        r_B_minus = 1.0 - r_B

        z_A = (np.log(r_A_plus) - np.log(r_A_minus)) / 2.0
        z_B = (np.log(r_B_plus) - np.log(r_B_minus)) / 2.0

        se = np.sqrt((1.0 / (n_A - 3)) + (1.0 / (n_B - 3)))
        z = (z_A - z_B) / se
        return z

    def z_to_p_value(self, z, two_sided=True):
        if two_sided:
            return scipy.stats.norm.sf(abs(z)) * 2
        else:
            return scipy.stats.norm.sf(abs(z))

    def build_miRNA_similarity_graph(self):
        self.G = copy.deepcopy(self.B.to_undirected())
        self.G.remove_nodes_from(nx.isolates(self.G))

        miRNAs_nodes = set(n for n, d in g.nodes(data=True) if d['bipartite'] == 0)
        targets_nodes = set(g) - miRNAs_nodes
        print 'mirnas', len(miRNAs_nodes)
        print 'targets', len(targets_nodes)
        print 'edges', len(g.edges())

        for m_i in miRNAs_nodes:
            for m_j in miRNAs_nodes:
                if (m_i is not m_j) and not (g.has_edge(m_i, m_j)):
                    common_neighbors = sorted(nx.common_neighbors(g, m_i, m_j))
                    if len(common_neighbors) > 0:
                        m_i_degree = self.G.degree(m_i)
                        m_j_degree = self.G.degree(m_j)

                        weight = 0.0
                        for gene in common_neighbors:
                            gene_degree = g.degree(gene)
                            weight += 1.0 / gene_degree * min(float(self.G.number_of_edges(m_i, gene)),
                                                              float(self.G.number_of_edges(m_j, gene)))

                        weight = min(weight / m_i_degree, weight / m_j_degree) + 1
                        if weight > 1.002:  #
                            self.G.add_edge(u=m_i, v=m_j, weight=weight)

        self.G.remove_nodes_from(targets_nodes)
        self.G.remove_nodes_from(nx.isolates(self.G))
        print len(set(self.G))
        print len(self.G.edges())

    def get_miRNA_community_assgn(self):
        partition = community.best_partition(self.G, weight='weight')
        # p_, nodes_community = zip(*sorted(partition.items()))
        return partition

    def build_miRNA_features(self, tags):
        dys_gene = OrderedDict()
        for tag in tags:
            dys_gene[tag] = np.unique([tup[1] for tup in self.B.edges(data=True) if tup[2]['tag'] == tag]).tolist()

        miRNAs_in_MTDN = np.unique([tup[0] for tup in self.B.edges()])
        print 'miRNAs_in_MTDN', len(miRNAs_in_MTDN)
        for tag, genes in dys_gene.iteritems():
            print 'genes dysregulated in ', tag, len(genes)

        genes_names = []
        for tag, genes in dys_gene.iteritems():
            for gene in genes:
                genes_names.append(gene + '/' + tag)

        self.miRNA_target_assn_matrix = pandas.DataFrame(columns=genes_names)

        for miRNA in miRNAs_in_MTDN:
            self.miRNA_target_assn_matrix.loc[miRNA] = 0

        for edge in self.B.edges(data=True):
            def f(x):
                dict = {}
                for tag, genes in dys_gene.iteritems():
                    if len(genes):
                        dict[tag] = 1.0 / len(genes)
                return dict[x]

            self.miRNA_target_assn_matrix.loc[edge[0]][edge[1] + '/' + edge[2]['tag']] = f(edge[2]['tag'])

    def run_miRNA_clustering(self, n_cluster=20, affinity='manhattan', linkage='complete'):

        mirna_cluster = AgglomerativeClustering(n_clusters=n_cluster, affinity=affinity, linkage=linkage).fit(
            self.miRNA_target_assn_matrix)
        self.miRNA_cluster_assgn = mirna_cluster.fit_predict(self.miRNA_target_assn_matrix)

        self.miRNA_clusters_int = []
        for cluster_idx in range(n_cluster):
            self.miRNA_clusters_int.append([])
            for mirna in [self.miRNA_target_assn_matrix.index[mirna_idx]
                          for mirna_idx, cluster_assg in enumerate(self.miRNA_cluster_assgn) if
                          cluster_assg == cluster_idx]:
                self.miRNA_clusters_int[cluster_idx].append(self.mirna_list.index(mirna))

        return np.bincount(mirna_cluster.labels_)

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
            if group_assg == -1 or group_assg == 0:
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
