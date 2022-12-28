#!/usr/bin/env python

from ase.calculators.vasp import Vasp
import ase.io
from ase.parallel import *
from assign_calc import base_calc
import numpy as np
import os
import argparse


def create_neb_directory_structure(atoms_is_sorted, atoms_fs_sorted, atoms_path_sorted):
    my_images = len(atoms_path_sorted)

    # Create directory structure for neb
    for i in np.arange(my_images+2):
        dir_name = str(i).zfill(2)
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        if dir_name == '00':
            ase.io.write(dir_name+'/'+'POSCAR', atoms_is_sorted,
                         vasp5=True, direct=True)
        elif dir_name == str(my_images+1).zfill(2):
            ase.io.write(dir_name+'/'+'POSCAR', atoms_fs_sorted,
                         vasp5=True, direct=True)
        else:
            atoms = atoms_path_sorted[i-1]
            ase.io.write(dir_name+'/'+'POSCAR', atoms, vasp5=True, direct=True)


def copy_contcar_to_poscar(my_images):
    current_dir = os.getcwd()
    os.system('mv INCAR INCAR.pre-neb')
    for i in np.arange(my_images):
        dir_name = str(i+1).zfill(2)
        # print("Copying CONTCAR to POSCAR for ", dir_name)
        os.chdir(dir_name)
        os.system(
            'cp CONTCAR CONTCAR.pre-neb; mv CONTCAR POSCAR; mv OUTCAR OUTCAR.pre-neb; mv OSZICAR OSZICAR.pre-neb')
        os.chdir(current_dir)


def presort(atoms):
    symbols = []
    symbolcount = {}
    for m, atom in enumerate(atoms):
        symbol = atom.symbol
        if symbol not in symbols:
            symbols.append(symbol)
            symbolcount[symbol] = 1
        else:
            symbolcount[symbol] += 1
    sort = []
    for symbol in symbols:
        for m, atom in enumerate(atoms):
            if atom.symbol == symbol:
                sort.append(m)
    resort = list(range(len(sort)))
    for n in range(len(resort)):
        resort[sort[n]] = n
    # atoms_sorted = atoms[sort]
    file = open('pre-sort.dat', 'w')  # same as ase-sort.dat
    for n in range(len(sort)):
        file.write('%5i %5i \n' % (sort[n], resort[n]))
    return sort


def run_initial_neb(atoms, my_images, my_nsw):
    """Run pre-neb for a given my_nsw steps"""

    atoms.calc.set(images=my_images, ibrion=3, iopt=0, lclimb=False,
                   spring=-5.00, potim=0.05, nsw=my_nsw)
    try:
        _ = atoms.get_potential_energy()  # Run vasp here
    except:
        IOError("Intended, as no CONTCAR is present. Continue to go to the actual neb.")


def run_actual_neb(atoms, my_images, my_nsw):
    """Run actual neb for a given my_nsw steps"""

    atoms.calc.set(images=my_images, ibrion=3, iopt=7, lclimb=True,
                   spring=-5.00, potim=0.0, nsw=my_nsw)

    try:
        _ = atoms.get_potential_energy()  # Run vasp here
    except:
        IOError("Intended, as no CONTCAR is present. Done calculation.")


def main(args):
    atoms_is = ase.io.read('is.traj')
    atoms_is.pbc = True
    atoms_fs = ase.io.read('fs.traj')
    atoms_fs.pbc = True
    guess_path = ase.io.read("guess_path.traj", index=":")

    if not os.path.exists("./all_neb_images.traj"):
        atoms_path = guess_path

    else:
        atoms_path = ase.io.read(
            "./all_neb_images.traj", index=f"-{len(guess_path)}:")

    # Sorting using ase.calculators.vasp.create_input
    sort = presort(atoms_is)
    atoms_is_sorted = atoms_is[sort]
    atoms_fs_sorted = atoms_fs[sort]
    atoms_path_sorted = []

    for atoms in atoms_path:
        atoms.pbc = True
        atoms_sorted = atoms[sort]
        atoms_path_sorted.append(atoms_sorted)

    # create dir structure if not present and writes POSCAR files
    create_neb_directory_structure(
        atoms_is_sorted, atoms_fs_sorted, atoms_path_sorted)

    my_images = len(atoms_path)

    atoms = atoms_is.copy()
    atoms.calc = base_calc()
    atoms.calc.set(encut=args.encut, ediffg=args.ediffg,
                   isif=2, xc="PBE", ivdw=12)

    if not os.path.exists("./all_neb_images.traj"):
        run_initial_neb(atoms, my_images, my_nsw=args.nsw_pre_neb)
        # Copy contcar to poscar for actual NEB
        copy_contcar_to_poscar(my_images)

    run_actual_neb(atoms, my_images, my_nsw=args.nsw_actual_neb)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=f"Runs VTST NEB calculation")
    parser.add_argument("--calc_dir", type=str, default=".",
                        help="path of NEB calculation directory (default is \".\" i.e. current working directory)", metavar="")
    parser.add_argument("--encut", type=float, default=400,
                        help="ENCUT to run the NEB (default is 400 eV)", metavar="")
    parser.add_argument("--ediffg", type=float, default=-0.05,
                        help="ediffg to run NEB (default is -0.05 eV/Ang)", metavar="")
    parser.add_argument("--nsw_pre_neb", type=int, default=5,
                        help="nsw for pre-NEB calculation (default is 5)", metavar="")
    parser.add_argument("--nsw_actual_neb", type=int, default=100,
                        help="nsw for actual-NEB calculation (default is 100)", metavar="")

    args = parser.parse_args()
    if args.calc_dir[-1] == "/":
        args.calc_dir = args.calc_dir[:-1]

    if not os.path.exists(args.calc_dir):
        raise f"Given path of calculation directory does not exist."
    else:
        start_dir = os.getcwd()
        os.chdir(args.calc_dir)
        main(args)
        os.chdir(start_dir)
