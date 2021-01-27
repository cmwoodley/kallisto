# tests/test_strucreader.py

import os

from atom import Atom
import reader.strucreader as ksr


# turbomole coord files can be read
def test_a_user_can_read_a_TMfile():
    f = open("tmpcoord", "w")
    s = os.linesep
    f.write("$coord" + s)
    f.write("  1.87167924 -0.101043656  0.1596818582  c" + s)
    f.write("  4.43543289 -0.101043656  0.1596818582  f" + s)
    f.write("  1.20847986 -1.386321988  1.6312493924  h" + s)
    f.write("  1.20847986 -0.732816897 -1.6891694984  h" + s)
    f.write("  1.20846096  1.816007916  0.5369845779  h" + s)
    f.write("$end")
    f.close()
    fileObject = open("tmpcoord", "r+")
    atoms = ksr.read(fileObject)
    os.remove("tmpcoord")
    got = type(atoms[0])
    want = Atom
    assert got is want
