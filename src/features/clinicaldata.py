import os

import pandas as pd

from definitions import ROOT_DIR


class ClinicalData:
    clinical_patient_colsname = ['bcr_patient_barcode', 'gender', 'race', 'histologic_diagnosis.1',
                                 'ajcc_pathologic_tumor_stage'
                                 ]

    pathologic_stage_map = {'Stage IA': 'Stage I', 'Stage IB': 'Stage I',
                            'Stage IIA': 'Stage II', 'Stage IIB': 'Stage II',
                            'Stage IIIA': 'Stage III', 'Stage IIIB': 'Stage III'}

    biospecimen_sample_colsname = ['bcr_sample_barcode', 'sample_type']

    clinical_drug_colsname = ['bcr_patient_barcode', 'pharmaceutical_therapy_drug_name', 'pharmaceutical_therapy_type',
                              'treatment_best_response']

    def __init__(self, cancer_type, folder_path):
        self.cancer_type = cancer_type

        # Import patients
        self.patient = pd.read_table(os.path.join(folder_path, "nationwidechildrens.org_clinical_patient.txt"),
                                     sep="\t",
                                     skiprows=[1, 2],
                                     na_values=["[Not Available]", "[Not Applicable]"],
                                     usecols=ClinicalData.clinical_patient_colsname
                                     )
        self.patient.index = self.patient["bcr_patient_barcode"]

        self.patient.replace({'ajcc_pathologic_tumor_stage': ClinicalData.pathologic_stage_map}, inplace=True)


        # # Import biospecimen samples (not all samples included in dataset)
        # self.biospecimen = pd.read_table(os.path.join(folder_path, "genome.wustl.edu_biospecimen_sample.txt"),
        #                                  sep="\t",
        #                                  skiprows=[1, ],
        #                                  na_values="[Not Available]",
        #                                  usecols=ClinicalData.biospecimen_sample_colsname
        #                                  )
        # self.biospecimen.index = self.biospecimen["bcr_sample_barcode"]

        # Import clinical drug
        self.drugs = pd.read_table(os.path.join(folder_path, "nationwidechildrens.org_clinical_drug.txt"),
                                   sep="\t",
                                   skiprows=[1, 2],
                                   na_values=["[Not Available]", "[Unknown]", "[Not Applicable]"],
                                   usecols=ClinicalData.clinical_drug_colsname
                                   )
        self.drugs.index = self.drugs["bcr_patient_barcode"]

        # Save index's
        self.patient_barcodes = self.patient["bcr_patient_barcode"].tolist()
        # self.sample_barcodes = self.biospecimen["bcr_sample_barcode"].tolist()
        # self.drug_barcodes = self.biospecimen["bcr_sample_barcode"].tolist()

    def get_patient_barcodes(self):
        return self.patient_barcodes

    # def get_sample_barcodes(self):
    #     return self.sample_barcodes


if __name__ == '__main__':
    # table = pd.read_table(ROOT_DIR+"/data/tcga-assembler/LUAD/clinical/nationwidechildrens.org_clinical_patient_luad.txt", sep="\t")
    folder_path = "/data/tcga-assembler/LUAD/clinical/"
    luad_clinical = ClinicalData(cancer_type="LUAD", folder_path=ROOT_DIR + folder_path)
