import os
import xml.etree.ElementTree as ET
import csv

"""
We extract relevant data fields from clinical xml files and combine into one csv file in /data/processed/clinical/
"""

clinical_src_dir = os.getcwd() + "/assn-mirna-luad/data/interim/clinical/"
clinical_tgt_dir = os.getcwd() + "/assn-mirna-luad/data/processed/clinical/"

luad = "http://tcga.nci/bcr/xml/clinical/luad/2.7"
shared = "http://tcga.nci/bcr/xml/shared/2.7"
clin_shared = "http://tcga.nci/bcr/xml/clinical/shared/2.7"
admin = "http://tcga.nci/bcr/xml/administration/2.7"
shared_stage = "http://tcga.nci/bcr/xml/clinical/shared/stage/2.7"
lung_shared = "http://tcga.nci/bcr/xml/clinical/shared/lung/2.7"

csv_file = open(clinical_tgt_dir + 'clinical.csv', 'w')
csvwriter = csv.writer(csv_file)
csvwriter.writerow(('patient_barcode', 'diagnosis', 'histological_type', 'pathologic_stage', 'pathologic_T',
                    'pathologic_N', 'pathologic_M'))

for file in os.listdir(clinical_src_dir):
    tree = ET.parse(clinical_src_dir + file)
    root = tree.getroot()

    for tcga_bcr in root.iter('{%s}tcga_bcr' % luad):
        patient = tcga_bcr.find('{%s}patient' % luad)
        histological_type = patient.find('{%s}histological_type' % shared).text
        patient_barcode = patient.find('{%s}bcr_patient_barcode' % shared).text
        diagnosis = patient.find('{%s}diagnosis' % lung_shared).text

        stage_event = patient.find('{%s}stage_event' % shared_stage)
        pathologic_stage = stage_event.find('{%s}pathologic_stage' % shared_stage).text

        tnm_categories = stage_event.find('{%s}tnm_categories' % shared_stage)
        pathologic_categories = tnm_categories.find('{%s}pathologic_categories' % shared_stage)
        pathologic_T = pathologic_categories.find('{%s}pathologic_T' % shared_stage).text
        pathologic_N = pathologic_categories.find('{%s}pathologic_N' % shared_stage).text
        pathologic_M = pathologic_categories.find('{%s}pathologic_M' % shared_stage).text

        if patient_barcode:
            row = patient_barcode, diagnosis, histological_type, pathologic_stage, pathologic_T, pathologic_N, pathologic_M
            # print row
            csvwriter.writerow(row)
        else:
            print patient_barcode, diagnosis, histological_type, pathologic_stage, pathologic_T, pathologic_N, pathologic_M