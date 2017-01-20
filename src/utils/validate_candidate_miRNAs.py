import os
from definitions import ROOT_DIR

VALIDATION_FILE_PATH = os.path.join(ROOT_DIR, 'data/external/TarBase_Experiment_Valid_miRNA-Targets.csv')

def percent_candidate_in_validated(candidate_miRNAs, validated_miRNAs):
    """
    Verify a model by percentages of known miRNAs found in experimentally
     validated databases (e.g. miRTarBase)

    :param candidate_miRNAs: a list of selected miRNAs code names
    :param validated_miRNAs: a list of validated miRNAs code names
    :return: a percentage
    """
    return float(len(set(candidate_miRNAs) & set(validated_miRNAs))) / len(validated_miRNAs)

def intersection_miRNA(candidate_miRNAs, validated_miRNAs):
    return set(candidate_miRNAs) & set(validated_miRNAs)