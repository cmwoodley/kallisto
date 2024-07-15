# tests/test_eeq.py
import numpy as np
from tests.store import ch_radical
from tests.store import acetylene
from tests.store import hydrogenCyanide

from kallisto.methods import get_cm5_corrections
from kallisto.data import CM5_ATOMIC_PARAMETERS
from kallisto.data import ATOMIC_RCOV

def test_eeq():
    charge = 0
    mol = ch_radical()
    eeq = mol.get_eeq(charge)
    assert np.isclose(eeq[0], -0.17166856)
    assert np.isclose(eeq[1], 0.17166856)


def test_eeq_cation():
    charge = 1
    mol = ch_radical()
    eeq = mol.get_eeq(charge)
    assert np.isclose(eeq[0], 0.59769359)
    assert np.isclose(eeq[1], 0.40230641)


def test_eeq_anion():
    charge = -1
    mol = ch_radical()
    eeq = mol.get_eeq(charge)
    assert np.isclose(eeq[0], -0.94103071)
    assert np.isclose(eeq[1], -0.05896929)


def test_eeq_cm5():
    charge = 0
    mol = ch_radical()
    cm5_charges = False
    eeq = mol.get_eeq(charge, cm5_charges)
    assert np.isclose(eeq[0], -0.17166856127379196)

    cm5_charges = True
    eeq = mol.get_eeq(charge, cm5_charges)
    assert np.isclose(eeq[0], -0.21562044659991492)

    # CM5 Correction Values as published in Table S4 J. Chem. Theory Comput. 2012, 8, 2, 527–541
    # DOI: 10.1021/ct200866d
    
    mol = acetylene()
    want = np.array([-0.053, -0.053, 0.053, 0.053])
    cm5_correction = mol.get_eeq(charge, True) - mol.get_eeq(charge, False)
    assert np.allclose(cm5_correction, want, atol=0.002)

    mol = hydrogenCyanide()
    want = np.array([0.067, -0.126, 0.059])
    cm5_correction = mol.get_eeq(charge, True) - mol.get_eeq(charge, False)
    assert np.allclose(cm5_correction, want, atol=0.002)


def test_cm5_atomic_parameters_access():
    # Test accessing CM5_ATOMIC_PARAMETERS by symbol
    i_sym = "H"
    j_sym = "Si"
    Dkk = CM5_ATOMIC_PARAMETERS[i_sym] - CM5_ATOMIC_PARAMETERS[j_sym]
    assert Dkk == 0.0846

    Dkk = CM5_ATOMIC_PARAMETERS[j_sym] - CM5_ATOMIC_PARAMETERS[i_sym]
    assert Dkk == -0.0846   

    i_sym = "H"
    j_sym = "N"
    assert (i_sym, j_sym) in CM5_ATOMIC_PARAMETERS.keys()
    Dkk = CM5_ATOMIC_PARAMETERS[(i_sym, j_sym)]
    assert Dkk == 0.1747

def test_atomic_rcov_access():
    # Test accessing ATOMIC_RCOV
    i_sym = "I"
    assert ATOMIC_RCOV[i_sym] == 2.57

    i_sym = "Be"
    assert ATOMIC_RCOV[i_sym] == 1.871


def test_get_cm5_correction():
    # Test CM5 corrections for CH radical
    mol = ch_radical()
    at = mol.get_atomic_numbers()
    coords = mol.get_positions()
    CM5_correction = get_cm5_corrections(at, coords)
    want = np.array([-0.04395187,  0.04395187])
    assert np.allclose(CM5_correction, want)

    # Test CM5 parameters for single hydrogen atom
    at = np.array([1])
    coords = np.array([[0,0,0]])
    CM5_correction = get_cm5_corrections(at, coords)
    assert np.isclose(CM5_correction[0], 0)

    # Test CM5 parameters for molecular nitrogen
    at = np.array([7, 7])
    # Example coords
    coords = np.array([[0., 0., 0.],
                       [1.0, 0., 0.]])
    CM5_correction = get_cm5_corrections(at, coords)
    want = np.array([0,0])
    assert np.allclose(CM5_correction, want)

    # Test CM5 parameters for distant atoms
    at = np.array([6, 1])
    # Example coords
    coords = np.array([[0., 0., 0.],
                       [20.0, 0., 0.]])
    CM5_correction = get_cm5_corrections(at, coords)
    want = np.array([0,0])
    assert np.allclose(CM5_correction, want)
