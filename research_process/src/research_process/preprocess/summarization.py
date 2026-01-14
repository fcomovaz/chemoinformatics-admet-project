from research_process.configs.configs_path import *
from research_process.utils.pytdc_endpoints import *

from research_process.configs.configs_log import log_flow, log_task
import logging

import pandas as pd

logger = logging.getLogger(__name__)


def perform_post_checks(
    filename: Path = DATASETS_DATA_DIR / "combined_admet.csv",
) -> None:
    """
    ## Perform post checks on the combined_admet.csv file.

    * Check if there are any duplicates in the file.
    * If there are duplicates, log a warning and write the duplicates to a file named "repeated_mols.csv" in the Processed Data directory.
    * Count the number of times each molecule appears in each process and write the counts to a file named "repeated_mols_summary.txt" in the Processed Data directory.

    :param filename: Path to the combined_admet.csv file

    :return: None
    """
    log_flow("Performing post checks - Duplicates")

    perform_repeat_reporting = False  # control variable, don't change

    df = pd.read_csv(filename)
    total_mols = df["Drug"]
    unique_mols = df["Drug"].unique()
    no_of_reps = len(total_mols) - len(unique_mols)

    if no_of_reps == 0:
        logger.info("No duplicates found in 'combined_admet.csv'")
    else:
        logger.warning(f"{no_of_reps} duplicates found in 'combined_admet.csv'")
        perform_repeat_reporting = True

    if perform_repeat_reporting:
        total_mols_rep = total_mols.duplicated(keep=False)
        total_mols_rep_len = len(total_mols_rep.values)

        set_mols_rep = set(total_mols[total_mols_rep])
        set_mols_rep_len = len(set_mols_rep)

        if set_mols_rep_len != total_mols_rep_len:
            logger.warning(f"Molecules ({set_mols_rep_len}) appears >= twice")
        logger.warning(f"Saving the {set_mols_rep_len} repeated molecules.")

        with open(PROCESSED_DATA_DIR / "repeated_mols.csv", "w") as f:
            f.write("SMILES\n")
            for mol in set_mols_rep:
                f.write(f"{mol}\n")
        logger.info(f"Molecules saved to {PROCESSED_DATA_DIR / 'repeated_mols.csv'}")

        # some numbers for the repetitions
        with open(PROCESSED_DATA_DIR / "repeated_mols_summary.txt", "w") as f:
            for name, file in admet_files.items():
                df = pd.read_csv(DATASETS_DATA_DIR / file)
                df_drug = df["Drug"]
                df_drug_len = df["Drug"].unique()
                total_repeats = len(df_drug) - len(df_drug_len)
                if total_repeats > 0:
                    logger.warning(f"Process {name:12} has {total_repeats:3} repeats")
                f.write(f"Process {name:12} has {total_repeats:3} repeats\n")

    logger.info("Post checks completed.")
