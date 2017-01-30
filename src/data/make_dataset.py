# -*- coding: utf-8 -*-
import os
import pandas
from definitions import ROOT_DIR
import numpy as np


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

        self.mirna_tumor.dropna(inplace=True)
        self.mirna_normal.dropna(inplace=True)

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

        self.gene_tumor.dropna(inplace=True)
        self.gene_normal.dropna(inplace=True)

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

        print self.mirna_tumor.shape, self.mirna_tumor.columns
        print self.mirna_normal.shape, self.mirna_normal.columns
        print self.gene_tumor.shape, self.gene_tumor.columns
        print self.gene_normal.shape, self.gene_normal.columns

    def make_dataset(self, dataset="miRNA", normal_tumor='both', pathologic_stages=[], normal_matched=True,
                     mirna_gene_matched=True, label_mapping=None):
        patients = pandas.concat([self.gene_tumor[['patient_barcode', 'pathologic_stage']],
                                  self.gene_normal[['patient_barcode', 'pathologic_stage']],
                                  self.mirna_tumor[['patient_barcode', 'pathologic_stage']],
                                  self.mirna_normal[['patient_barcode', 'pathologic_stage']]]).drop_duplicates()
        print 'All patients', patients.shape

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

        if mirna_gene_matched and normal_matched:
            if dataset is "miRNA":
                patients_B = pandas.merge(self.mirna_tumor[['patient_barcode']],
                                          self.mirna_tumor[['patient_barcode']],
                                          on='patient_barcode')
            elif dataset is "gene":
                patients_B = pandas.merge(self.gene_tumor[['patient_barcode']],
                                          self.gene_tumor[['patient_barcode']],
                                          on='patient_barcode')

            patients = pandas.merge(patients, patients_B, on='patient_barcode')
        elif mirna_gene_matched and not normal_matched:
            print "Not Implemented"
            return None

        print "normal_matched", normal_matched, 'mirna_gene_matched', mirna_gene_matched, dataset, normal_tumor, patients.shape

        if dataset is 'miRNA':
            if normal_tumor is 'both':
                return self.dataFrame_to_matrix(pandas.concat([self.mirna_tumor, self.mirna_normal]), patients,
                                                pathologic_stages, label_mapping)
            elif normal_tumor is 'normal':
                return self.dataFrame_to_matrix(self.mirna_normal, patients,
                                                pathologic_stages, label_mapping)
            elif normal_tumor is 'tumor':
                return self.dataFrame_to_matrix(self.mirna_tumor, patients,
                                                pathologic_stages, label_mapping)
        elif dataset is 'gene':
            if normal_tumor is 'both':
                return self.dataFrame_to_matrix(pandas.concat([self.gene_tumor, self.gene_normal]), patients,
                                                pathologic_stages, label_mapping)
            elif normal_tumor is 'normal':
                return self.dataFrame_to_matrix(self.gene_normal, patients,
                                                pathologic_stages, label_mapping)
            elif normal_tumor is 'tumor':
                return self.dataFrame_to_matrix(self.gene_tumor, patients,
                                                pathologic_stages, label_mapping)

    def dataFrame_to_matrix(self, data, patients, pathologic_stages, label_mapping):
        print 'data.shape', data.shape
        df = data[data['patient_barcode'].isin(patients['patient_barcode'])]
        print 'df.shape', df.shape

        if pathologic_stages:
            df = df[df['pathologic_stage'].isin(pathologic_stages)]

        if label_mapping:
            df = df['pathologic_stage'].replace(label_mapping)

        X = df.drop(['patient_barcode', 'pathologic_stage'], axis=1)
        y = df['pathologic_stage']

        return X, y

    def get_miRNA_list(self):
        return self.mirna_list

    def get_gene_list(self):
        return self.gene_symbols


if __name__ == '__main__':
    tgca_luad = TCGA_LUAD()
    X, y = tgca_luad.make_dataset(dataset='miRNA', normal_tumor='both', normal_matched=False, mirna_gene_matched=False,
                                  pathologic_stages=[])
    print "X", X.shape
    print "y", y.shape
