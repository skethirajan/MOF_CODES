#!/usr/bin/env python

import sys
import ase.atoms
import ase.io
from ase.parallel import *
from assign_calc import isif_2_calc


def run_isif_2(atoms: ase.atoms.Atoms) -> None:
    calc = isif_2_calc()
    atoms.calc = calc
    _ = atoms.get_potential_energy()


if __name__ == '__main__':
    file = sys.argv[1]
    if ".traj" in file:
        atoms_1 = ase.io.read(filename=file, index=-1, format="traj")
    elif ".xyz" in file:
        atoms_1 = ase.io.read(filename=file, index=-1, format="xyz")
    elif ".xml" in file:
        atoms_1 = ase.io.read(filename=file, index=-1, format="vasp-xml")
    else:
        raise f"Check the format of the file, {file} to read."
    run_isif_2(atoms=atoms_1)
