from rdkit.Chem import Descriptors

TOTAL_DESCRIPTORS = len(Descriptors.descList)

ALL_DESCRIPTORS = [desc for desc in Descriptors.descList]

AVG_TIME = 150 / 8000  # seconds
