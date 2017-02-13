import os
from operator import itemgetter

from definitions import ROOT_DIR

VALIDATION_FILE_PATH = os.path.join(ROOT_DIR, 'data/external/TarBase_Experiment_Valid_miRNA-Targets.csv')


def recall_rate(candidate_miRNAs, validated_miRNAs):
    """
    Measures recall rate. Percent of candidate miRNA's selected, in validated miRNA's

    :param candidate_miRNAs: a list of selected miRNAs code names
    :param validated_miRNAs: a list of validated miRNAs code names
    :return: a percentage
    """
    return float(len(intersection_miRNA(candidate_miRNAs, validated_miRNAs))) / len(validated_miRNAs)


def precision_rate(candidate_miRNAs, validated_miRNAs):
    """
    Measures precision rate. Percent of validated miRNA's selected, in candidate miRNA's

    :param candidate_miRNAs: a list of selected miRNAs code names
    :param validated_miRNAs: a list of validated miRNAs code names
    :return: a percentage
    """
    return float(len(intersection_miRNA(candidate_miRNAs, validated_miRNAs))) / len(candidate_miRNAs)

def intersection_miRNA(candidate_miRNAs, validated_miRNAs):
    return set(candidate_miRNAs) & set(validated_miRNAs)


def get_miRNA_names(indices, mirna_list, miR_name=False):
    if miR_name:
        result = list(itemgetter(*indices)(mirna_list))
        for i in range(len(result)):
            result[i] = result[i].replace('hsa-', '')
        return result
    return itemgetter(*indices)(mirna_list)
