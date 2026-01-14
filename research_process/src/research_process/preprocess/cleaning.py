from research_process.configs.configs_path import *
from research_process.utils.pytdc_endpoints import *

from research_process.configs.configs_log import log_flow, log_task
import logging

import pandas as pd

from rdkit import Chem

logger = logging.getLogger(__name__)


def removing_duplicates() -> None:
    """
    Removes duplicates from the combined_admet.csv file.

    * Reads the repeated_mols.csv file and the combined_admet.csv file.
    * If the repeated_mols.csv file is not found, an empty dataframe is loaded.
    * If the combined_admet.csv file is not found, the function ends.
    * Renames the "Drug" column to "SMILES".
    * Removes the duplicates from the combined_admet.csv file based on the repeated_mols.csv file.
    * Saves the new dataframe to combined_admet_no_dups.csv.

    :param: None

    :return: None
    """
    log_flow("Removing duplicates")

    logger.info("Loading repeated_mols.csv file")
    try:
        mols = pd.read_csv(PROCESSED_DATA_DIR / "repeated_mols.csv")
    except FileNotFoundError:
        logger.error("File 'repeated_mols.csv' not found.")
        mols = pd.DataFrame(columns=["SMILES"])
        logger.warning("No duplicates found. Loading empty dataframe.")

    logger.info("Loading combined_admet.csv file")
    try:
        df = pd.read_csv(DATASETS_DATA_DIR / "combined_admet.csv")
    except FileNotFoundError:
        logger.error("File 'combined_admet.csv' not found.")
        return

    try:
        df = df.rename(columns={"Drug": "SMILES"})
    except Exception as e:
        logger.error(f"Error renaming columns: {e}")

    df = df[~df["SMILES"].isin(mols["SMILES"])]

    df.to_csv(DATASETS_DATA_DIR / "combined_admet_no_dups.csv", index=False)

    logger.info("Duplicates removed")


def removing_non_valid_smiles() -> None:
    """
    Removes non valid SMILES from the combined_admet_no_dups.csv file.

    * Reads the combined_admet_no_dups.csv file.
    * If the combined_admet_no_dups.csv file is not found, the function ends.
    * Removes the non valid SMILES from the combined_admet_no_dups.csv file.
    * Saves the new dataframe to combined_admet_no_invalids.csv.

    :param: None

    :return: None
    """
    log_flow("Removing non valid SMILES - rdkit check")

    logger.info("Loading combined_admet_no_dups.csv file")
    try:
        df = pd.read_csv(DATASETS_DATA_DIR / "combined_admet_no_dups.csv")
    except FileNotFoundError:
        logger.error("File 'combined_admet_no_dups.csv' not found.")
        return

    logger.info(f"Dataset size before removal {df.shape[0]:,}")

    logger.info("Removing non valid SMILES")
    TOTAL_INVALID = 0
    invalid_smiles = []
    mols = df["SMILES"].tolist()
    for m in mols:
        _m = Chem.MolFromSmiles(m)
        if _m is None:
            logger.error(f"Invalid SMILES: {m}")
            TOTAL_INVALID += 1
            invalid_smiles.append(m)
            df = df[df["SMILES"] != m]

    logger.info(f"Total invalid SMILES: {TOTAL_INVALID}")
    logger.info(f"Dataset size after  removal {df.shape[0]:,}")

    with open(PROCESSED_DATA_DIR / "invalid_smiles.txt", "w") as f:
        f.write("\n".join(invalid_smiles))

    df.to_csv(DATASETS_DATA_DIR / "combined_admet_no_invalids.csv", index=False)

    logger.info("Non valid SMILES removed")
