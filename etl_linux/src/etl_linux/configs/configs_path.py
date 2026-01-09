from pathlib import Path


def find_repo_root(marker="data"):
    p = Path(__file__).resolve()
    for parent in p.parents:
        if (parent / marker).exists():
            return parent
    raise RuntimeError("Repo root not found")


REPO_ROOT = find_repo_root("data")  # find the root, not the data/ folder

RESULTS_DIR = REPO_ROOT / "results"

DATA_DIR = REPO_ROOT / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
ZIPPED_DATA_DIR = DATA_DIR / "zipped"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
