import os

import numpy as np
import pandas as pd

from definitions import ROOT_DIR


class GenomicData(object):
    def __init__(self, cancer_type, file_path, columns="GeneSymbol|TCGA", log2_transform=True):
        self.cancer_type = cancer_type

        self.data = self.preprocess_expression_table(pd.read_table(file_path, sep="\t"), columns)

        if log2_transform:
            self.data = self.data.applymap(self.log2_transform)

        # Save samples and features for this omics data
        self.samples = self.data.index
        self.features = self.data.columns.tolist()
        # self.features.remove("bcr_sample_barcode")

    def preprocess_expression_table(self, df, columns):
        """
        Download
        :param df:
        :param columns:
        :return:
        """
        table = df

        # Filter columns
        table = table.filter(regex=columns)

        # Cut TCGA column names to sample barcode
        table.rename(columns=lambda x: x[:16] if ("TCGA" in x) else x, inplace=True)

        # Drop duplicate columns names (Gene symbols with same name)
        _, i = np.unique(table.columns, return_index=True)
        table = table.iloc[:, i]

        # Drop NA GeneSymbol rows
        table.dropna(axis=0, inplace=True)

        # Remove entries with unknown Gene Symbol
        table = table[table.GeneSymbol != '?']

        # Transpose dataframe to patient rows and GeneSymbol columns
        table.index = table.GeneSymbol
        table.drop(['GeneSymbol'], axis=1, inplace=True)
        table = table.T

        # Add column for patients barcode
        # table['bcr_sample_barcode'] = table.index

        # Drop duplicate columns names (Gene symbols with same name)
        _, i = np.unique(table.columns, return_index=True)
        table = table.iloc[:, i]

        return table

    def log2_transform(self, x):
        return np.log2(x+1)


    def get_genes_list(self):
        return self.features

    def get_samples_list(self):
        return self.samples


class LncRNAExpression(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "TCGA-rnaexpr.tsv")
        super(LncRNAExpression, self).__init__(cancer_type, file_path)

    def preprocess_expression_table(self, df, columns):
        lncrna_exp = df

        lncrna_names = pd.read_table(
            os.path.join(ROOT_DIR, "data/tcga-assembler/LUAD/lncrna/HGNC_RNA_long_non-coding.txt"),
            delimiter="\t")
        lncrna_dict = pd.Series(lncrna_names.symbol.values, index=lncrna_names.ensembl_gene_id).to_dict()

        # Replacing ENSG Gene ID to the lncRNA symbol name
        lncrna_exp['Gene_ID'] = lncrna_exp['Gene_ID'].str.replace("[.].*", "")
        lncrna_exp.replace({"Gene_ID": lncrna_dict}, inplace=True)

        # Drop NA gene rows
        lncrna_exp.dropna(axis=0, inplace=True)

        # Transpose matrix to patients rows and genes columns
        lncrna_exp.index = lncrna_exp['Gene_ID']
        lncrna_exp = lncrna_exp.T.iloc[1:, :]

        # Change index string to bcr_sample_barcode standard
        def change_patient_barcode(s):
            if "Normal" in s:
                return s[s.find('TCGA'):] + "-11A"
            elif "Tumor" in s:
                return s[s.find('TCGA'):] + "-01A"
            else:
                return s

        lncrna_exp.index = lncrna_exp.index.map(change_patient_barcode)

        return lncrna_exp


class GeneExpression(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "geneExp.txt")
        super(GeneExpression, self).__init__(cancer_type, file_path)


class SNP(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "somaticMutation_geneLevel.txt")
        super(SNP, self).__init__(cancer_type, file_path)


class miRNAExpression(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "miRNAExp__RPM.txt")
        super(miRNAExpression, self).__init__(cancer_type, file_path)


class CopyNumberVariation(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "copyNumber.txt")
        super(CopyNumberVariation, self).__init__(cancer_type, file_path)


class DNAMethylation(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "")
        super(DNAMethylation, self).__init__(cancer_type, file_path)


class ProteinExpression(GenomicData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "protein_RPPA.txt")
        super(ProteinExpression, self).__init__(cancer_type, file_path)


if __name__ == '__main__':
    # table = pd.read_table(ROOT_DIR+"/data/tcga-assembler/LUAD/clinical/nationwidechildrens.org_clinical_patient_luad.txt", sep="\t")
    folder_path = "/data/tcga-assembler/LUAD/lncrna/"
    lncRNA_expression = LncRNAExpression(cancer_type="LUAD", folder_path=ROOT_DIR + folder_path)
