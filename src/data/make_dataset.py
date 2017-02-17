# -*- coding: utf-8 -*-
import os

import numpy as np
import pandas
from definitions import ROOT_DIR
from sklearn.feature_selection import SelectFdr, f_classif


class TCGA_LUAD:
    def __init__(self):
        # miRNA
        mirna_tumor_df = pandas.read_csv(os.path.join(ROOT_DIR, "data/processed/miRNA/tumor_miRNA.csv"))
        mirna_normal_df = pandas.read_csv(os.path.join(ROOT_DIR, "data/processed/miRNA/normal_miRNA.csv"))
        self.clinical_df = pandas.read_csv(os.path.join(ROOT_DIR, "data/processed/clinical/clinical.csv"))

        self.mirna_normal = pandas.merge(self.clinical_df[['patient_barcode', 'pathologic_stage']], mirna_normal_df,
                                         on='patient_barcode')
        self.mirna_normal['pathologic_stage'] = 'normal'
        self.mirna_tumor = pandas.merge(self.clinical_df[['patient_barcode', 'pathologic_stage']], mirna_tumor_df,
                                        on='patient_barcode')

        self.mirna_tumor.dropna(axis=0, inplace=True)
        self.mirna_normal.dropna(axis=0, inplace=True)

        pathologic_stage_map = {'Stage IA': 'Stage I', 'Stage IB': 'Stage I',
                                'Stage IIA': 'Stage II', 'Stage IIB': 'Stage II',
                                'Stage IIIA': 'Stage III', 'Stage IIIB': 'Stage III'}
        self.mirna_tumor.replace({'pathologic_stage': pathologic_stage_map}, inplace=True)

        self.mirna_list = list(self.mirna_tumor.columns)[2:]

        # Gene Expression
        gene_tumor_df = pandas.read_table(
            os.path.join(ROOT_DIR, 'data/processed/gene_expression/tumor/READ__illuminahiseq_rnaseqv2__GeneExp.txt'),
            header=0, delimiter='\t')
        gene_normal_df = pandas.read_table(
            os.path.join(ROOT_DIR, 'data/processed/gene_expression/normal/READ__illuminahiseq_rnaseqv2__GeneExp.txt'),
            header=0, delimiter='\t')

        gene_tumor_df.rename(columns=lambda x: x[:12], inplace=True)
        gene_normal_df.rename(columns=lambda x: x[:12], inplace=True)

        # Remove entries with unknown Gene Symbol
        gene_tumor_df = gene_tumor_df[gene_tumor_df.GeneSymbol != '?']
        gene_normal_df = gene_normal_df[gene_normal_df.GeneSymbol != '?']

        # Get list of all gene_symbols
        self.gene_symbols = list(gene_tumor_df['GeneSymbol'])

        # Get list of tumor and normal patient_barcode
        gene_exp_tumor_patient_barcodes = list(gene_tumor_df.columns)[2:]
        gene_exp_normal_patient_barcodes = list(gene_normal_df.columns)[2:]

        # Drop EntrezID column
        self.gene_tumor = gene_tumor_df.drop(['EntrezID', 'GeneSymbol'], axis=1)
        self.gene_normal = gene_normal_df.drop(['EntrezID', 'GeneSymbol'], axis=1)

        # Reshaping data frame to have columns for GeneSymbols, and rows of patients
        self.gene_tumor = self.gene_tumor.T
        self.gene_normal = self.gene_normal.T
        self.gene_tumor.columns = self.gene_symbols
        self.gene_normal.columns = self.gene_symbols

        # Add column for patients barcode
        self.gene_tumor['patient_barcode'] = self.gene_tumor.index
        self.gene_normal['patient_barcode'] = self.gene_normal.index

        self.gene_tumor.dropna(axis=0, inplace=True)
        self.gene_normal.dropna(axis=0, inplace=True)

        self.gene_normal = pandas.merge(self.clinical_df[['patient_barcode', 'pathologic_stage']], self.gene_normal,
                                        on='patient_barcode')
        self.gene_normal['pathologic_stage'] = 'normal'
        self.gene_tumor = pandas.merge(self.clinical_df[['patient_barcode', 'pathologic_stage']], self.gene_tumor,
                                       on='patient_barcode')

        self.gene_tumor.replace({'pathologic_stage': pathologic_stage_map}, inplace=True)

        # Drop duplicate columns names (Gene symbols with same name)
        _, i = np.unique(self.gene_tumor.columns, return_index=True)
        self.gene_tumor = self.gene_tumor.iloc[:, i]
        self.gene_normal = self.gene_normal.iloc[:, i]

        print 'mirna_tumor', self.mirna_tumor.shape
        print 'mirna_normal', self.mirna_normal.shape
        print 'gene_tumor', self.gene_tumor.shape
        print 'gene_normal', self.gene_normal.shape

    def gene_univariate_feature_selection(self, alpha=0.01):
        gene_normal_X, gene_normal_Y = self.make_dataset(dataset='gene', normal_tumor='normal',
                                                         normal_matched=True, mirna_gene_matched=True)
        gene_tumor_X, gene_tumor_Y = self.make_dataset(dataset='gene', normal_tumor='tumor', normal_matched=True,
                                                       mirna_gene_matched=True)

        gene_exp_filter = SelectFdr(f_classif, alpha=alpha)
        gen_exp_new = gene_exp_filter.fit_transform(X=pandas.concat([gene_normal_X, gene_tumor_X]),
                                                    y=pandas.concat([gene_normal_Y, gene_tumor_Y]))

        self.gene_symbols = np.asanyarray(self.gene_symbols)[gene_exp_filter.get_support(indices=True)].tolist()
        self.gene_tumor = self.gene_tumor[self.gene_symbols + ['patient_barcode', 'pathologic_stage']]
        self.gene_normal = self.gene_normal[self.gene_symbols + ['patient_barcode', 'pathologic_stage']]


    def make_dataset(self, dataset="miRNA", normal_tumor='both', pathologic_stages=[], normal_matched=True,
                     mirna_gene_matched=True, label_mapping=None, zero_mean=False, normalize=False):

        # Find patients with both tumor and normal samples
        if normal_matched:
            if dataset is "miRNA":
                patients = pandas.merge(self.mirna_tumor[['patient_barcode']],
                                        self.mirna_normal[['patient_barcode']],
                                        on='patient_barcode')
            elif dataset is "gene":
                patients = pandas.merge(self.gene_tumor[['patient_barcode']],
                                        self.gene_normal[['patient_barcode']],
                                        on='patient_barcode')
        elif not normal_matched:
            if dataset is "miRNA":
                patients = pandas.concat([self.mirna_tumor[['patient_barcode']],
                                          self.mirna_normal[['patient_barcode']]]).drop_duplicates()
            elif dataset is "gene":
                patients = pandas.concat([self.gene_tumor[['patient_barcode']],
                                          self.gene_normal[['patient_barcode']]]).drop_duplicates()

        # Find patients with matching miRNA and gene samples
        if mirna_gene_matched and normal_matched:
            patients_normal = pandas.merge(self.mirna_normal[['patient_barcode']],
                                           self.gene_normal[['patient_barcode']],
                                           on='patient_barcode')
            patients_tumor = pandas.merge(self.mirna_tumor[['patient_barcode']],
                                          self.gene_tumor[['patient_barcode']],
                                          on='patient_barcode')
            patients = pandas.merge(patients_normal, patients_tumor, on='patient_barcode')

        elif mirna_gene_matched and not normal_matched:
            if normal_tumor is "normal":
                patients = pandas.merge(self.mirna_normal[['patient_barcode']],
                                        self.gene_normal[['patient_barcode']],
                                        on='patient_barcode')
            elif normal_tumor is "tumor":
                patients = pandas.merge(self.mirna_tumor[['patient_barcode']],
                                        self.gene_tumor[['patient_barcode']],
                                        on='patient_barcode')

        # Return dataset, and perform pathologic stage relabling
        if dataset is 'miRNA':
            if normal_tumor is 'both':
                return self.dataFrame_to_matrix(pandas.concat([self.mirna_tumor, self.mirna_normal]), patients,
                                                pathologic_stages, label_mapping, zero_mean, normalize)
            elif normal_tumor is 'normal':
                return self.dataFrame_to_matrix(self.mirna_normal, patients,
                                                pathologic_stages, label_mapping, zero_mean, normalize)
            elif normal_tumor is 'tumor':
                return self.dataFrame_to_matrix(self.mirna_tumor, patients,
                                                pathologic_stages, label_mapping, zero_mean, normalize)
        elif dataset is 'gene':
            if normal_tumor is 'both':
                return self.dataFrame_to_matrix(pandas.concat([self.gene_tumor, self.gene_normal]), patients,
                                                pathologic_stages, label_mapping, zero_mean, normalize)
            elif normal_tumor is 'normal':
                return self.dataFrame_to_matrix(self.gene_normal, patients,
                                                pathologic_stages, label_mapping, zero_mean, normalize)
            elif normal_tumor is 'tumor':
                return self.dataFrame_to_matrix(self.gene_tumor, patients,
                                                pathologic_stages, label_mapping, zero_mean, normalize)

    def dataFrame_to_matrix(self, data, patients, pathologic_stages, label_mapping, zero_mean=False, normalize=False):
        df = data[data['patient_barcode'].isin(patients['patient_barcode'])]
        if pathologic_stages:
            df = df[df['pathologic_stage'].isin(pathologic_stages)]

        if label_mapping:
            df['pathologic_stage'] = df['pathologic_stage'].replace(label_mapping)

        X = df.drop(['patient_barcode', 'pathologic_stage'], axis=1)
        y = df['pathologic_stage']

        if normalize:
            for col in X.columns:
                X[col] = (X[col] - X[col].mean()) / X[col].std(0)
        elif zero_mean:
            for col in X.columns:
                X[col] = X[col] - X[col].mean()

        return X, y

    def get_miRNA_list(self):
        return self.mirna_list

    def get_gene_list(self):
        return self.gene_symbols


if __name__ == '__main__':
    tgca_luad = TCGA_LUAD()
    # pathologic_stage_map = {'Stage I': 1, 'Stage II': 1}
    X, y = tgca_luad.make_dataset(dataset='miRNA', normal_tumor='tumor', normal_matched=False, mirna_gene_matched=True,
                                  pathologic_stages=[])
    print "X", X.shape
    print "y", y.shape

    print tgca_luad.make_dataset(dataset='gene', normal_tumor='tumor', normal_matched=False, mirna_gene_matched=True)[
        0].shape
