import pandas as pd
import numpy as np

from definitions import ROOT_DIR
from src.features.clinicaldata import ClinicalData
from src.features.genomicdata import GeneExpression, SNP, DNAMethylation, miRNAExpression, CopyNumberVariation, \
    ProteinExpression, LncRNAExpression
# from src.features.slide_image import WholeSlideImages


class MultiOmicsData:
    pathologic_stage_map = {'Stage IA': 'Stage I', 'Stage IB': 'Stage I',
                            'Stage IIA': 'Stage II', 'Stage IIB': 'Stage II',
                            'Stage IIIA': 'Stage III', 'Stage IIIB': 'Stage III'}

    def __init__(self, cancer_type, folder_path, modalities=["WSI", "GE", "SNP", "CNV", "DNA", "MIR", "PRO"]):
        """
        Load all multi-omics TCGA data from a given folder_path with the following folder structure:

            folder_path/
                clinical/
                    genome.wustl.edu_biospecimen_sample.txt
                    nationwidechildrens.org_clinical_patient.txt
                gene_exp/
                    geneExp.txt
                mirna/
                    miRNAExp__RPM.txt
                cnv/
                    copyNumber.txt
                protein_rppa/
                    protein_RPPA.txt
                somatic/
                    somaticMutation_geneLevel.txt

        :param cancer_type: TCGA cancer code name
        :param folder_path: relative directory path to the folder containing multi-omics data downloaded from TCGA-assembler
        """
        self.cancer_type = cancer_type
        self.modalities = modalities

        # LOADING DATA FROM FILES
        self.multi_omics_data = {}
        self.clinical = ClinicalData(cancer_type, folder_path + "clinical/")
        self.multi_omics_data["PATIENTS"] = self.clinical.patient
        # self.multi_omics_data["BIOSPECIMENS"] = self.clinical.biospecimen
        self.multi_omics_data["DRUGS"] = self.clinical.drugs

        if ("WSI" in modalities):
            self.WSI = WholeSlideImages(cancer_type, folder_path)
            self.multi_omics_data["WSI"] = self.WSI
        if ("GE" in modalities):
            self.GE = GeneExpression(cancer_type, folder_path + "gene_exp/", )
            self.multi_omics_data["GE"] = self.GE.data
        if ("SNP" in modalities):
            self.SNP = SNP(cancer_type, folder_path + "somatic/")
            self.multi_omics_data["SNP"] = self.SNP.data
        if ("MIR" in modalities):
            self.MIR = miRNAExpression(cancer_type, folder_path + "mirna/")
            self.multi_omics_data["MIR"] = self.MIR.data
        if ("LNC" in modalities):
            self.LNC = LncRNAExpression(cancer_type, folder_path + "lncrna/")
            self.multi_omics_data["LNC"] = self.LNC.data
        if ("DNA" in modalities):
            self.DNA = DNAMethylation(cancer_type, folder_path + "dna/")
            self.multi_omics_data["DNA"] = self.DNA.data
        if ("CNV" in modalities):
            self.CNV = CopyNumberVariation(cancer_type, folder_path + "cnv/")
            self.multi_omics_data["CNV"] = self.CNV.data
        if ("PRO" in modalities):
            self.PRO = ProteinExpression(cancer_type, folder_path + "protein_rppa/")
            self.multi_omics_data["PRO"] = self.PRO.data

        self.print_sample_sizes()

    def match_samples(self, modalities):
        """
        Return the index of bcr_sample_barcodes of the intersection of samples from all modalities

        :param modalities: An array of modalities
        :return: An pandas Index list
        """
        matched_samples = self.multi_omics_data[modalities[0]].index.copy()

        for modality in modalities:
            matched_samples = matched_samples.join(self.multi_omics_data[modality].index, how="inner")

        return matched_samples

    def load_data(self, multi_omics, target=['ajcc_pathologic_tumor_stage', 'histologic_diagnosis.1'], pathologic_stages=[], histological_types=[]):
        """
        Load and return the multi-omics dataset (classification)
        :param multi_omics: A list of the data modalities to load. Default "all" to select all modalities
        """
        if multi_omics == 'all' or multi_omics == None:
            modalities = self.modalities
        elif multi_omics:
            modalities = multi_omics
        else:
            raise Exception("Need to specify which multi-omics to retrieve")

        matched_samples = self.match_samples(modalities)

        # Build targets clinical data
        y = self.get_patients_clinical(matched_samples)


        # Filter samples
        y = y[y['ajcc_pathologic_tumor_stage'] != "[Discrepancy]"]
        y.loc[y.index.str.contains("-11A"),
              'ajcc_pathologic_tumor_stage'] = "Normal"  # Change stage label of normal samples to "Normal"

        # Select only samples with certain cancer stage or subtype
        if pathologic_stages:
            y = y[y['ajcc_pathologic_tumor_stage'].isin(pathologic_stages)]
        if histological_types:
            y = y[y['histologic_diagnosis.1'].isin(histological_types)]

        # Filter y target column labels
        y = y.filter(target)
        y.dropna(axis=0, inplace=True)

        matched_samples = y.index

        # Build data matrix for each modality, indexed by matched_samples
        X_multiomics = {}
        for modality in modalities:
            X_multiomics[modality] = self.multi_omics_data[modality].loc[matched_samples]

        return X_multiomics, y

    def get_patients_clinical(self, matched_samples, clinical_data=['PATIENTS']):
        """
        Fetch patient's clinical data given a list of samples barcodes

        :param matched_samples:
        :param clinical_data:
        :return:
        """
        target = pd.DataFrame(index=matched_samples)

        # Make a separate column for patients barcode from samples barcode
        target["patient_barcode"] = target.index.str[:-4]
        print("modalities matched sample size:", target.shape)

        # Merge patients clinical data with patient barcode as index
        for clinical in clinical_data:
            target = target.join(self.multi_omics_data[clinical],
                                 how="left", on="patient_barcode", rsuffix="_")

        # TODO if normal_matched:
        #     target =
        print("joined clinical data size:", target.shape)

        target.replace({'ajcc_pathologic_tumor_stage': MultiOmicsData.pathologic_stage_map}, inplace=True)

        return target  # Return only the columns specified


    def print_sample_sizes(self):
        for modality in self.multi_omics_data.keys():
            print(modality, self.multi_omics_data[modality].shape if hasattr(self.multi_omics_data[modality],
                                                                             'shape') else "Didn't import data")



if __name__ == '__main__':
    folder_path = "/home/jonny_admin/PycharmProjects/nuclei-segmentation/data/tcga-assembler/LUSC/"
    luad_data = MultiOmicsData(cancer_type="LUSC", folder_path=folder_path,
                               modalities=["GE", "MIR"])

    matched_samples = luad_data.match_samples(modalities=["GE", "MIR"])
    print("matched samples", matched_samples.shape, matched_samples)

    # X,  = luad_data.load_data(multi_omics= "all", target=['ajcc_pathologic_tumor_stage'])
    patients_clinical = luad_data.get_patients_clinical(matched_samples)
    # print(patients_clinical)
