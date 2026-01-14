from research_process.configs.configs_log import setup_logging
from research_process.preprocess.summarization import *

setup_logging("research_process")

if __name__ == "__main__":
    perform_post_checks()
