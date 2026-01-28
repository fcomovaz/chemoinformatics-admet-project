# pretuning.py
from research_process.configs.configs_log import log_flow, log_task
from research_process.configs.configs_path import PROCESSED_DATA_DIR
import logging

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

import bblean
import bblean.plotting as plotting
import bblean.analysis as analysis

import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def bit_performance() -> None:
    """
    Evaluate the density of fingerprints at different bit lengths.

    Compute the density of fingerprints for the molecules in molecules_smiles.csv
    at different bit lengths (2^n where n is between 0 and 14).
    The density is computed as the ratio of the number of bits set to the total
    number of bits in the fingerprint.

    The results are logged to the console.

    Returns:
        None
    """
    # -------------------------------------------------
    # Config
    # -------------------------------------------------
    BITS = [2**n for n in range(15)]

    try:
        df = pd.read_csv(PROCESSED_DATA_DIR / "molecules_smiles.csv")
        df = df.sample(frac=0.1, random_state=42).reset_index(drop=True)
        mols = df["SMILES"].tolist()
    except FileNotFoundError:
        logger.error("File 'molecules_smiles.csv' not found.")
        return
    try:
        mols = [Chem.MolFromSmiles(s) for s in mols if Chem.MolFromSmiles(s)]
    except Exception as e:
        logger.error(f"Error converting SMILES: {e}")
        return

    # -------------------------------------------------
    # FP density evaluation
    # -------------------------------------------------
    results = []

    for size in BITS:
        gen = GetMorganGenerator(radius=2, fpSize=size)

        # dens = []
        # for m in mols:
        #     fp = gen.GetFingerprint(m)
        #     dens = 100.0 * fp.GetNumOnBits() / size
        #     dens.append(dens)
        try:
            dens = [100.0 * gen.GetFingerprint(m).GetNumOnBits() / size for m in mols]
        except Exception as e:
            logger.error(f"Error computing density: {e}")
            return

        mean_density = np.mean(dens)
        std_density = np.std(dens)

        results.append((size, mean_density, std_density))

    for size, mean, std in results:
        logger.info(f"{size:>6} bits  →  {mean:6.2f}% ± {std:5.2f}%")

    return


def assign_cluster(clusters) -> np.ndarray:
    """
    Assigns each molecule to a cluster.

    Given a list of clusters (a list of lists of molecule indices),
    assigns each molecule to a cluster and returns a numpy array
    with the cluster assignments.

    Parameters
    ----------
    clusters : list of lists of int
        A list of clusters, where each cluster is a list of molecule indices.

    Returns
    -------
    numpy.ndarray
        A numpy array with the cluster assignments for each molecule.
    """
    log_flow("Assigning clusters")
    total_items = sum([len(c) for c in clusters])
    logger.info(f"Total items: {total_items}")

    cluster_list = np.full(total_items, np.nan)
    for i, cluster in enumerate(clusters):
        for mol_index in cluster:
            cluster_list[mol_index] = i

    return cluster_list


def preselect_bit_bitbirch() -> None:
    log_task("preselect_bit_bitbirch")

    try:
        df = pd.read_csv(PROCESSED_DATA_DIR / "molecules_smiles.csv")
        df = df.sample(frac=0.1, random_state=42).reset_index(drop=True)
        mols = df["SMILES"].tolist()
    except FileNotFoundError:
        logger.error("File 'molecules_smiles.csv' not found.")
        return
    # try:
    #     mols = [Chem.MolFromSmiles(s) for s in mols if Chem.MolFromSmiles(s)]
    # except Exception as e:
    #     logger.error(f"Error converting SMILES: {e}")
    # return

    fps = bblean.fps_from_smiles(mols, pack=True, n_features=4096, kind="ecfp4")

    tree = bblean.BitBirch(threshold=0.35, merge_criterion="tolerance-diameter")
    tree.fit(fps)
    clusters = tree.get_cluster_mol_ids()
    logger.info(f"BitBIRCH Clustering =========> {len(clusters)} clusters")

    inner_items = [len(c) for c in clusters]
    singletons = inner_items.count(1)
    logger.info(f"Singletons: {singletons/len(df)*100:.4f}%")
    ca = analysis.cluster_analysis(clusters, fps, mols)
    plotting.summary_plot(ca, title="Molecules - ECFP4 - BitBirch Clustering")
    plt.show()

    tree.set_merge(merge_criterion="tolerance-radius", tolerance=0.0)
    tree.refine_inplace(fps)
    clusters = tree.get_cluster_mol_ids()
    logger.info(f"Refined BitBIRCH Clustering => {len(clusters)} clusters")

    inner_items = [len(c) for c in clusters]
    singletons = inner_items.count(1)
    logger.info(f"Singletons: {singletons/len(df)*100:.4f}%")
    ca = analysis.cluster_analysis(clusters, fps, mols)
    plotting.summary_plot(ca, title="Molecules - ECFP4 - BitBirch Clustering")
    plt.show()
