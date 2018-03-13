import pandas as pd

from src.features.clinicaldata import ClinicalData
from src.features.genomicdata import GeneExpression, SNP, DNAMethylation, miRNAExpression, CopyNumberVariation, \
    ProteinExpression, LncRNAExpression


# from src.features.slide_image import WholeSlideImages


class MultiOmicsData:

    def __init__(self, cancer_type, folder_path, modalities=["WSI", "GE", "SNP", "CNV", "DNA", "MIR", "PRO"]):
        """
        Load all multi-omics TCGA data from a given folder_path with the following folder structure:

            folder_path/
                clinical/
                    genome.wustl.edu_biospecimen_sample.txt
                    nationwidechildrens.org_clinica l_patient.txt
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
        # TODO check that for single modalities, this fetch all patients
        matched_samples = self.multi_omics_data[modalities[0]].index.copy()

        for modality in modalities:
            matched_samples = matched_samples.join(self.multi_omics_data[modality].index, how="inner")

        return matched_samples

    def load_data(self, multi_omics, target=['ajcc_pathologic_tumor_stage'],
                  pathologic_stages=[], histological_types=[], predicted_subtypes=[], samples_barcode=None):
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
        if not (samples_barcode is None):
            matched_samples = samples_barcode

        # Build targets clinical data
        y = self.get_patients_clinical(matched_samples)
        # print(y)

        # Filter samples
        y = y[y['ajcc_pathologic_tumor_stage'] != "[Discrepancy]"]
        y.loc[y.index.str.contains("-11A"),
              'ajcc_pathologic_tumor_stage'] = "Normal"  # Change stage label of normal samples to "Normal"

        # Select only samples with certain cancer stage or subtype
        if pathologic_stages:
            y = y[y['ajcc_pathologic_tumor_stage'].isin(pathologic_stages)]
        if histological_types:
            y = y[y['histologic_diagnosis.1'].isin(histological_types)]
        if predicted_subtypes:
            y = y[y['predicted_subtype'].isin(predicted_subtypes)]
        # TODO if normal_matched:
        #     target =

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
        Fetch patient's clinical data for each given samples barcodes in matched_samples

        :param matched_samples:
        :param clinical_data:
        :return:
        """
        target = pd.DataFrame(index=matched_samples)

        # Make a separate column for patients barcode from samples barcode
        target["patient_barcode"] = target.index.str[:-4]

        no_samples = target.shape[0]

        # Merge patients clinical data with patient barcode as index
        for clinical in clinical_data:
            # target = pd.merge(target, self.multi_omics_data[clinical],
            #                      how="left", left_on="patient_barcode", right_on="patient_barcode")
            target = target.join(self.multi_omics_data[clinical], on="patient_barcode", how="left", rsuffix="_")

        # print(target)
        if target.shape[0] != no_samples:
            raise Exception("Clinical data merging has wrong number of samples")

        return target  # Return only the columns specified


    def print_sample_sizes(self):
        for modality in self.multi_omics_data.keys():
            print(modality, self.multi_omics_data[modality].shape if hasattr(self.multi_omics_data[modality],
                                                                             'shape') else "Didn't import data")

    def add_subtypes_to_patients_clinical(self, dictionary):
        self.multi_omics_data["PATIENTS"] = self.multi_omics_data["PATIENTS"].assign(
            predicted_subtype=self.multi_omics_data["PATIENTS"]["bcr_patient_barcode"].map(dictionary))


if __name__ == '__main__':
    folder_path = "/home/jonny_admin/PycharmProjects/nuclei-segmentation/data/tcga-assembler/LUSC/"
    luad_data = MultiOmicsData(cancer_type="LUSC", folder_path=folder_path,
                               modalities=["GE", "MIR"])

    matched_samples = luad_data.match_samples(modalities=["GE", "MIR"])
    print("matched samples", matched_samples.shape, matched_samples)

    # X,  = luad_data.load_data(multi_omics= "all", target=['ajcc_pathologic_tumor_stage'])
    patients_clinical = luad_data.get_patients_clinical(matched_samples)
    # print(patients_clinical)
