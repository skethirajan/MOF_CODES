#!/usr/bin/env python

import sys
import ase.atoms
import ase.io
from ase.parallel import *
from assign_calc import nvt_calc


def run_nvt(atoms: ase.atoms.Atoms, T_set: float) -> None:
    calc = nvt_calc()
    calc.set(tebeg=T_set)
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
    T = float(sys.argv[2])  # in Kelvin
    run_nvt(atoms=atoms_1, T_set=T)
