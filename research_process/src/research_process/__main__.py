from research_process.configs.configs_log import setup_logging

from research_process.preprocess.summarization import *
from research_process.preprocess.cleaning import *
from research_process.preprocess.enrichment import *
from research_process.clustering.pretuning import *

setup_logging("research_process")

if __name__ == "__main__":
    # summarization module
    # perform_post_checks()

    # cleaning module
    # removing_duplicates()
    # removing_non_valid_smiles()
    # check_descriptors_availability()
    # separate_final_molecules()
    # create_incomplete_endpoints()

    # enrichment module
    # perform_enrichment()

    # clustering module

    # bit_performance()
    preselect_bit_bitbirch()
