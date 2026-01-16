from research_process.configs.configs_path import *
from research_process.utils.pytdc_endpoints import *
from research_process.utils.rdkit_descriptors import *
from research_process.configs.configs_log import log_flow, log_task
import logging
import pandas as pd
import numpy as np
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


def check_descriptors_availability() -> None:
    """
    Checks if all molecules have valid descriptors.

    * Reads the combined_admet_no_invalids.csv file.
    * Computes all descriptors for each molecule.
    * Saves the molecules with valid descriptors to combined_admet_full_descrips.csv.
    * Saves the molecules with invalid descriptors to invalid_descriptors_molecules.txt.

    :param: None

    :return: None
    """
    log_flow("Checking descriptors availability")
    df = pd.read_csv(DATASETS_DATA_DIR / "combined_admet_no_invalids.csv")
    logger.info(f"Dataset size {df.shape[0]:,}")
    # df = df.sample(frac=0.01, random_state=42).reset_index(drop=True)
    mols = df["SMILES"].tolist()

    proc_time = AVG_TIME * len(mols)
    logger.warning(f"This process may take at least {proc_time:6.2f} seconds")
    logger.warning(f"This process may take at least {proc_time / 60:6.2f} minutes")

    descriptors = []
    for i, m in enumerate(mols):
        # print(f"Processed -> {i+1}/{len(mols)}")
        m = Chem.MolFromSmiles(m)
        if m is None:
            logger.error(f"Invalid SMILES: {m}")
            break

        try:
            descriptors.append({name: f(m) for name, f in ALL_DESCRIPTORS})
        except Exception:
            # Mol válida pero descriptor falló → NaN también
            descriptors.append({name: np.nan for name, _ in ALL_DESCRIPTORS})
            logger.error(f"Error computing descriptors: {m}")

    mols = df["SMILES"].reset_index(drop=True)
    df = pd.concat([mols, pd.DataFrame(descriptors)], axis=1)

    logger.info("Getting molecules with invalid descriptors")

    non_desc = df[df.isna().any(axis=1)]
    yes_desc = df[~df.isna().any(axis=1)]

    logger.info(f"Total molecules with invalid descriptors: {non_desc.shape[0]:,}")
    logger.info(f"Total molecules with valid descriptors: {yes_desc.shape[0]:,}")

    with open(PROCESSED_DATA_DIR / "invalid_descriptors_molecules.txt", "w") as f:
        f.write("\n".join(non_desc["SMILES"].tolist()))

    yes_desc.to_csv(DATASETS_DATA_DIR / "combined_admet_full_descrips.csv", index=False)

    logger.info("Invalid descriptors removed")


def separate_final_molecules() -> None:
    """
    Separates the final molecules from the combined_admet_full_descrips.csv file.

    * Reads the combined_admet_full_descrips.csv file.
    * Resets the index of the SMILES column.
    * Saves the SMILES column to final_mols.csv.

    :param: None

    :return: None
    """

    log_flow("Separating final molecules")

    df = pd.read_csv(DATASETS_DATA_DIR / "combined_admet_full_descrips.csv")
    logger.info(f"Dataset size {df.shape[0]:,}")

    mols = df["SMILES"]

    mols = mols.reset_index(drop=True)

    mols.to_csv(PROCESSED_DATA_DIR / "molecules_smiles.csv", index=False)

    logger.info("Final molecules separated")
