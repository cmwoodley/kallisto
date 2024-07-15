# optimisation/config.py

# directory configuration
DATA_PATH = "data"
MOLS_DIRNAME = "qm_charges"
OPTUNA_DIRNAME = "optuna_fitting"

# file configuration
ORIGINAL_DATA_FILENAME = "all_mols_def2TZVP_charges_cm5.sdf"
PKBHX_DATA_FILENAME = "pkbhx.txt"
OPTUNA_CM5_STUDY_FILENAME = "01-optuna_CM5_fitting_study.pbz2"
OPTUNA_CM5_PARAMETERS_FILENAME = "02-optuna_CM5_fitting_best_params.json"

# logic configuration
VERBOSE = True
OPTUNA_EARLY_STOPPING = 300
LOSS_FUNCTION = "mean_absolute_error"
FITTING_SMARTS = [
            "[NX2,nX3][OX1]",  # General N-O single bond
            "[N+][O-]", # Aliphatic N-oxide or Nitro, 
            "[n+][O-]", # Aromatic N-oxide
            "[S]=[O]", # S=O double bond
            "[P]=[O]",  # P=O double bond
        ]
CM5_PARAMETERS_UNMODIFIED= {
    ("H", "C"): 0.0502,
    ("H", "N"): 0.1747,
    ("H", "O"): 0.1671,
    ("C", "N"): 0.0556,
    ("C", "O"): 0.0234,
    ("N", "O"): -0.0346,
    "H": 0.0056,
    "He": -0.1543,
    "Li": 0.0,
    "Be": 0.0333,
    "B": -0.1030,
    "C": -0.0446,
    "N": -0.1072,
    "O": -0.0802,
    "F": -0.0629,
    "Ne": -0.1088,
    "Na": 0.0184,
    "Mg": 0.0,
    "Al": -0.0726,
    "Si": -0.0790,
    "P": -0.0756,
    "S": -0.0565,
    "Cl": -0.0444,
    "Ar": -0.0767,
    "K": 0.0130,
    "Ca": 0.0,
    "Zn": 0.0,
    "Ge": -0.0557,
    "As": -0.0533,
    "Se": -0.0399,
    "Br": -0.0313,
    "I": -0.0220,
    "Xe": -0.0381,
    "Cs": 0.0065,
    "Ba": 0,
    "La": 0,
    "Hf": 0,
    "Ta": 0,
    "W": 0,
    "Re": 0,
    "Os": 0,
    "Ir": 0,
    "Pt": 0,
    "Au": 0,
    "Hg": 0,
    "Tl": -0.0255,
    "Pb": -0.0277,
    "Bi": -0.0265,
    "Po": -0.0198,
    "At": -0.0155,
    "Rn": -0.0269,
    "Fr": 0.0046,
    "Ra": 0,
}
