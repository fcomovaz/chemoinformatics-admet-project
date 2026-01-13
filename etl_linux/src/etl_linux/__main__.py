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
    is_part_of = set(files).issubset(set(admet_files.values()))
    if is_part_of and len(files) > 0:
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


# # -----------------------------------------------------------
# # -----------------------------------------------------------
# log_flow("Collapsing over duplicates with different values")
# is_grouping_performed = False  # change to False to avoid grouping the variables
# if is_grouping_performed:
#     for admet_process in admet_processes:
#         log_task(f"Grouping {admet_process}")
#         try:
#             file_name = admet_files[admet_process]
#             df = pd.read_csv(DATASETS_DATA_DIR / file_name)
#         except Exception as e:
#             logger.error(f"Error reading file '{file_name}': {e}")

#         try:
#             # df = df.groupby("Drug").agg(lambda x: x.mode().iloc[0])
#             df = df.groupby("Drug", as_index=False).mean()
#             df.to_csv(DATASETS_DATA_DIR / file_name, index=False)
#         except Exception as e:
#             logger.error(f"Error grouping over duplicates: {e}")

#         logger.info(f"Grouping completed. File {file_name} processed.")
# else:
#     logger.info("Grouping not performed.")


# -----------------------------------------------------------
# -----------------------------------------------------------
log_flow("Creating the big dataset - Joint ADMET")
big_df = None
for idx, admet_process in enumerate(admet_processes):
    log_task(f"Processing {admet_process} ({idx+1}/{len(admet_processes)})")

    try:
        file_name = admet_files[admet_process]
        df = pd.read_csv(DATASETS_DATA_DIR / file_name)
    except Exception as e:
        logger.error(f"Error reading file '{file_name}': {e}")

    # prefix columns with process name to avoid duplicates
    df = df.rename(
        columns={c: f"{admet_process}_{c}" if c != "Drug" else c for c in df.columns}
    )

    if big_df is None:
        big_df = df
    else:
        big_df = pd.merge(big_df, df, on="Drug", how="outer")

    logger.info(f"Processing completed. File '{file_name}' processed.")


big_df.to_csv(DATASETS_DATA_DIR / "combined_admet.csv", index=False)
log_flow(f"Combined DataFrame saved to '{DATASETS_DATA_DIR}/combined_admet.csv'")
