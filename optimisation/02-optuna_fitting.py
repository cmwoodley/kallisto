"""NO, PO, SO parameter fitting using Optuna and sklearn."""

# optimisation/optuna_fitting.py
import config
from helpers import get_cm5_corrections_dict
from helpers import kallisto_molecule_from_rdkit_molecule
from helpers import MolsFromSpartanFiles
from helpers import get_SMARTS_matches
from helpers import run_optimisation
from helpers import dump_parameters_to_json
import os
import numpy as np

"""
Notes:
The parameter optimisation was carried out using optuna==3.6.1
(https://pypi.org/project/optuna/)

[I 2024-07-12 12:25:48,111] Trial 526 finished with value: 0.11406869094467403 and parameters: {'NO': 0.26429166589238134, 'SO': -0.05372604653179257, 'PO': -0.12288464204139021}. Best is trial 224 with value: 0.1136626516981626.
Early stopping exceeded: No new best scores             after 300 iterations
The best parameters were {'NO': 0.2734730210431726, 'SO': -0.06094900180684683, 'PO': -0.14257153352901053}
"""

# Import load files for input
data_path = os.path.abspath(os.path.join(os.getcwd(), "..", config.DATA_PATH))
optuna_path = os.path.join(data_path, config.OPTUNA_DIRNAME)
study_filepath = os.path.join(optuna_path, config.OPTUNA_CM5_STUDY_FILENAME)
params_filepath = os.path.join(optuna_path, config.OPTUNA_CM5_PARAMETERS_FILENAME)
original_dataset_filepath = os.path.join(
    data_path, config.MOLS_DIRNAME, config.ORIGINAL_DATA_FILENAME
)

# Load in molecules with MP2/6-31G* CHELPG charges
mols = MolsFromSpartanFiles(original_dataset_filepath)

# Identify molecules affected by NO, PO and SO CM5 pairwise corrections
smarts_matches = get_SMARTS_matches(config.FITTING_SMARTS, mols)

print(len(smarts_matches))

# Pre-calculate uncorrected kallisto charges and extract QM charges
input_data = dict()
for idx in smarts_matches:
    mol = mols[idx]
    mol_data = dict()
    mol_data["rdkit_mol"] = mol
    kmol = kallisto_molecule_from_rdkit_molecule(mol)
    mol_data["kallisto_mol"] = kmol
    mol_data["qm_charges"] = np.array(
        # [float(atom.GetProp("_TriposPartialCharge")) for atom in mol.GetAtoms()]
        [float(x) for x in mol.GetProp("chelpg_charges").split(",")]
    )
    mol_data["kallisto_charges"] = kmol.get_eeq(0)
    mol_data["at"] = kmol.get_atomic_numbers()
    mol_data["coords"] = kmol.get_positions()

    input_data[idx] = mol_data


def objective(trial):

    # Define the equation parameters to optimize
    NO = trial.suggest_float("NO", -0.5, 0.5)
    SO = trial.suggest_float("SO", -0.5, 0.5)
    PO = trial.suggest_float("PO", -0.5, 0.5)

    # Set up updated CM5 corrections dictionary
    cm5_dict = config.CM5_PARAMETERS_UNMODIFIED
    cm5_dict[("N", "O")] = NO
    cm5_dict[("O", "P")] = PO
    cm5_dict[("O", "S")] = SO

    # Calculate per-atom error
    errors = []
    for idx in input_data:
        at = input_data[idx]["at"]
        coords = input_data[idx]["coords"]
        kallisto_charges = input_data[idx]["kallisto_charges"]
        qm_charges = input_data[idx]["qm_charges"]
        CM5_correction = get_cm5_corrections_dict(
            at, coords, CM5_ATOMIC_PARAMETERS=cm5_dict
        )
        cm5_charges = kallisto_charges + CM5_correction

        errors.append(cm5_charges - qm_charges)

    # metric to minimise
    if config.LOSS_FUNCTION == "mean_absolute_error":
        loss = np.mean(np.abs(np.concatenate(errors)))
    elif config.LOSS_FUNCTION == "mean_squared_error":
        loss = np.mean(np.concatenate(errors) ** 2)
    elif config.LOSS_FUNCTION == "root_mean_square_error":
        loss = np.sqrt(np.mean(np.concatenate(errors) ** 2))
    else:
        raise ValueError("Please define a valid loss function")
    return loss


# Run optimisation
run_optimisation(objective, study_filepath, verbose=True)

# Dump best parameters to json
dump_parameters_to_json(study_filepath, params_filepath)
