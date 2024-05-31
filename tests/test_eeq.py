# tests/test_eeq.py
import numpy as np
from tests.store import ch_radical
from tests.store import acetylene
from tests.store import hydrogenCyanide


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

    mol = acetylene()
    want = np.array([-0.053, -0.053, 0.053, 0.053])
    cm5_correction = mol.get_eeq(charge, True) - mol.get_eeq(charge, False)
    assert np.allclose(cm5_correction, want, atol=0.002)

    mol = hydrogenCyanide()
    want = np.array([0.067, -0.126, 0.059])
    cm5_correction = mol.get_eeq(charge, True) - mol.get_eeq(charge, False)
    assert np.allclose(cm5_correction, want, atol=0.002)
