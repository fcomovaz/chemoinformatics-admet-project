from research_process.configs.configs_path import *
from research_process.utils.pytdc_endpoints import *

from research_process.configs.configs_log import log_flow, log_task
import logging

import pandas as pd

logger = logging.getLogger(__name__)


def perform_enrichment() -> None:
    log_flow("Performing enrichment - Creating Descriptors Dataset")

    df = pd.read_csv(DATASETS_DATA_DIR / "combined_admet_full_descrips.csv")
    logger.info(f"Dataset size {df.shape[0]:,}")

    df.to_csv(PROCESSED_DATA_DIR / "molecules_with_descriptors.csv", index=False)

    logger.info("Enrichment separation performed")
