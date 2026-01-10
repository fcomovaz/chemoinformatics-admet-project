# etl_linux/__main__.py
from etl_linux.configs.configs_path import RAW_DATA_DIR, DATASETS_DATA_DIR
from etl_linux.utils.pytdc_endpoints import *
from tdc.single_pred import ADME, Tox  # database from https://tdcommons.ai/single-pred

from etl_linux.configs.configs_log import *
import logging

import pandas as pd
import os


setup_logging("etl")  # Logging configuration
logger = logging.getLogger(__name__)  # start logging


# -----------------------------------------------------------
# -----------------------------------------------------------
is_caching_necessary = False  # DON'T CHANGE
try:
    files = os.listdir(RAW_DATA_DIR)
    total_files = len(files) - 2  # -2 for the .gitkeep/dataset folder
    logger.info(f"Folder {RAW_DATA_DIR} contains {total_files} files.")
    if total_files == TOTAL_ENDPOINTS:
        logger.info("Caching is not necessary.")
    else:
        logger.warning("Caching is necessary.")
        is_caching_necessary = True
except Exception as e:
    logger.error(f"Error getting {RAW_DATA_DIR} files: {e}")


# -----------------------------------------------------------
# -----------------------------------------------------------
if is_caching_necessary:
    try:
        for admet_process in admet_processes:
            log_flow(f"Caching files for {admet_process} at {RAW_DATA_DIR}")
            for admet_property in admet_processes[admet_process]:
                log_task(f"Processing {admet_property}")
                if admet_process == "Toxicity":
                    tox = Tox(name=admet_property, path=RAW_DATA_DIR)
                else:
                    adme = ADME(name=admet_property, path=RAW_DATA_DIR)
    except Exception as e:
        logger.error(f"Error caching: {e}")


# -----------------------------------------------------------
# -----------------------------------------------------------
try:
    if os.path.exists(DATASETS_DATA_DIR):
        logger.info(f"Folder {DATASETS_DATA_DIR} detected")
    else:
        logger.warning(f"No such folder {DATASETS_DATA_DIR}.")
        os.mkdir(DATASETS_DATA_DIR)
        logger.info(f"Folder {DATASETS_DATA_DIR} created")
except Exception as e:
    logger.error(f"Error creating {DATASETS_DATA_DIR} folder: {e}")


# -----------------------------------------------------------
# -----------------------------------------------------------
is_acquisition_necessary = False  # DON'T CHANGE
try:
    # get the list from datasets
    files = os.listdir(DATASETS_DATA_DIR)
    if set(files).issubset(set(admet_files.values())):
        logger.info("Acquisition is not necessary.")
    else:
        logger.warning("Acquisition is necessary.")
        is_acquisition_necessary = True
except Exception as e:
    logger.error(f"Error getting {DATASETS_DATA_DIR} files: {e}")


# -----------------------------------------------------------
# -----------------------------------------------------------
if is_acquisition_necessary:
    try:
        for admet_process in admet_processes:
            log_flow(f"Processing files for {admet_process} from {RAW_DATA_DIR}")
            try:
                file_name = admet_files[admet_process]
                cols = [prprty for prprty in admet_processes[admet_process]]
                cols.append("Drug")
                # logger.info(f"Creating columns: {cols}")
            except Exception as e:
                logger.error(f"Error creating columns: {e}")

            df = None
            for admet_property in admet_processes[admet_process]:
                log_task(f"Processing {admet_property}")
                try:
                    if admet_process == "Toxicity":
                        tox = Tox(name=admet_property, path=RAW_DATA_DIR)
                        data = tox.get_data()
                    else:
                        adme = ADME(name=admet_property, path=RAW_DATA_DIR)
                        data = adme.get_data()
                except Exception as e:
                    logger.error(f"Error getting data: {e}")

                try:
                    # convert data to pandas
                    data_df = pd.DataFrame(data)
                except Exception as e:
                    logger.error(f"Error converting data to pandas: {e}")

                try:
                    # rename Y by admet_property
                    data_df.rename(columns={"Y": admet_property}, inplace=True)
                except Exception as e:
                    logger.error(f"Error renaming columns: {e}")

                try:
                    # just two columns should survive ["Drug", admet_property]
                    survivor_cols = ["Drug", admet_property]
                    data_df = data_df[survivor_cols]
                except Exception as e:
                    logger.error(f"Error preserving columns: {e}")

                try:
                    # add admet_property column
                    if df is None:
                        df = data_df  # first time
                    else:
                        # df = pd.merge(df, data_df, on='Drug', how='inner') # no nan
                        df = pd.merge(df, data_df, on="Drug", how="outer")  # allow nan
                except Exception as e:
                    logger.error(f"Error merging dataframes {admet_property}: {e}")

            # drop duplicates (e.g. A=5, A=5)
            df = df.drop_duplicates()  # but not double data (e.g. A=5, A=2)
            df.to_csv(DATASETS_DATA_DIR / file_name, index=False)

            logger.info(f"Processing completed. File '{file_name}' created.")
    except Exception as e:
        logger.error(f"Error processing files: {e}")
