# Single-instance Prediction ./Datasets
absorption = [
    "Caco2_Wang",
    "HIA_Hou",
    "Pgp_Broccatelli",
    "Bioavailability_Ma",
    "Lipophilicity_AstraZeneca",
    "Solubility_AqSolDB",
    "HydrationFreeEnergy_FreeSolv",
]
distribution = [
    "BBB_Martins",
    "PPBR_AZ",
    "VDss_Lombardo",
]
metabolism = [
    "CYP2C19_Veith",
    "CYP2D6_Veith",
    "CYP3A4_Veith",
    "CYP_1A2_Veith",
    "CYP2C9_Veith",
    "CYP2C9_Substrate_CarbonMangels",
    "CYP2D6_Substrate_CarbonMangels",
    "CYP3A4_Substrate_CarbonMangels",
]
execretion = [
    "Half_Life_Obach",
    "Clearance_Hepatocyte_AZ",
    "Clearance_Microsome_AZ",
]
toxicity = [
    "hERG",
    "AMES",
    "DILI",
    "Skin Reaction",
    "LD50_Zhu",
    "Carcinogens_Lagunin",
    # "ToxCast", # not working for some reason
    # "Tox21",  # not working for some reason
    "ClinTox",  # maybe not useful
]

# Dictionary for each process
admet_processes = {
    "Absorption": absorption,
    "Distribution": distribution,
    "Metabolism": metabolism,
    "Excretion": execretion,
    "Toxicity": toxicity,
}
