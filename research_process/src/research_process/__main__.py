from research_process.configs.configs_log import setup_logging

from research_process.preprocess.summarization import *
from research_process.preprocess.cleaning import *
from research_process.preprocess.enrichment import *

setup_logging("research_process")

if __name__ == "__main__":
    # perform_post_checks()

    # removing_duplicates()
    # removing_non_valid_smiles()
    # check_descriptors_availability()
    separate_final_molecules()

    # perform_enrichment()
