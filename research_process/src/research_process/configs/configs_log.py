# configs_log.py
import logging as lg
import functools
from datetime import datetime
from pathlib import Path

FLOW_LINE = "#" * 60
TASK_LINE = "-" * 60


def _build_log_path(name: str) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H")
    log_dir = Path(__file__).resolve().parents[2] / "logs"
    log_dir.mkdir(exist_ok=True)
    return log_dir / f"{name}_{timestamp}.log"


def setup_logging(filename="logs.log", level=lg.INFO):
    log_path = _build_log_path(filename)

    root_logger = lg.getLogger()  # global logger
    root_logger.setLevel(level)  # set global log level - script level

    # If it already has handlers, don't add more
    if root_logger.handlers:
        return root_logger

    # Format
    # formatter = lg.Formatter("[%(asctime)s][%(levelname)s] - %(message)s")
    # formatter = lg.Formatter("[%(asctime)s][%(levelname)-8s] - %(message)s")
    formatter = lg.Formatter("[%(levelname).1s][%(asctime)s][%(name)s] - %(message)s")
    # formatter = lg.Formatter("[%(asctime)s][%(levelname)s][%(name)s] - %(message)s")

    # File handler
    fh = lg.FileHandler(log_path, mode="a", encoding="utf-8-sig")
    fh.setFormatter(formatter)

    # Console handler
    ch = lg.StreamHandler()
    ch.setFormatter(formatter)

    # Add handlers to the logger
    root_logger.addHandler(fh)
    root_logger.addHandler(ch)

    return root_logger


def flow(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        logger = lg.getLogger(func.__module__)

        logger.info("#" * 60)

        try:
            return func(*args, **kwargs)
        finally:
            logger.info("#" * 60)

    return wrapper


def log_flow(msg: str, logger=None):
    logger = logger or lg.getLogger()
    logger.info(FLOW_LINE)
    logger.info(msg)
    logger.info(FLOW_LINE)


def log_task(msg: str, logger=None):
    logger = logger or lg.getLogger()
    logger.info(TASK_LINE)
    logger.info(msg)
    logger.info(TASK_LINE)


# minimum working example
# from {folder}.configs.configs_log import *
# import logging

# setup_logging("etl")  # Logging configuration
# logger = logging.getLogger(__name__)  # start logging
