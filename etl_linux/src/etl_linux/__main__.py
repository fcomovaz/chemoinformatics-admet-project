# etl_linux/__main__.py
from etl_linux.configs.configs_log import *
import logging

from etl_linux.utils.pytdc_endpoints import *

setup_logging("etl")  # Logging configuration
logger = logging.getLogger(__name__)  # start logging

log_flow("This is a flow example")
logger.info("This is a step 1")
logger.info("This is a step 2")
log_task("This is a task example")
logger.info("This is a step 3")
logger.error("This is an error message")
logger.critical("HOLA")
logger.warning(f"{absorption}")

from etl_linux.configs.configs_path import *

# save a file into that dir
logger.info(f"WE ARE IN {REPO_ROOT}")
with open(ZIPPED_DATA_DIR / "test.txt", "w") as f:
    f.write("Hello World")
