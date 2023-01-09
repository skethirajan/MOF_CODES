#!/usr/bin/env python

import ase.atoms
import ase.io
from ase.parallel import *
from assign_calc import base_calc
import argparse
import os
import re


def run_isif_2(atoms: ase.atoms.Atoms, kwargs) -> None:
    calc = base_calc()
    calc.set(isif=2, nsw=kwargs.nsw, nelm=kwargs.nelm)
    atoms.calc = calc
    _ = atoms.get_potential_energy()


def get_atoms(file, index):
    if ".traj" in file:
        atoms = ase.io.read(filename=file, index=index, format="traj")
    elif ".xyz" in file:
        atoms = ase.io.read(filename=file, index=index, format="xyz")
    elif ".xml" in file:
        atoms = ase.io.read(filename=file, index=index, format="vasp-xml")
    else:
        raise f"Check the format of the file, {file} to read."
    return atoms


def read_command_line_arguments():
    parser = argparse.ArgumentParser(description="Runs VASP isif=2 relaxation")
    parser.add_argument("file", type=str, help="file path of the structure")
    parser.add_argument(
        "--index", type=int, default=-1, help="index of the structure to consider in the given file (default is last index)", metavar="")
    parser.add_argument("--nsw", type=int, default=500,
                        help="number of ionic steps (default is 500)", metavar="")
    parser.add_argument("--nelm", type=int, default=100,
                        help="max. electronic self-consistency steps (default is 100)", metavar="")
    parser.add_argument("--new_lattice", type=float,
                        help="set the cell with new lattice value", metavar="")
    parser.add_argument("--restart", action="store_true",
                        help="upon restart, structure from last index in vasprun.xml will be used")
    args = parser.parse_args()
    return args


def print_verbose(args):
    with open("./job.out", "a") as fi:
        fi.write(f"Ionic steps remaining: {args.nsw}\n")
        fi.write(
            f"Using structure at index={args.index} from file: {args.file}\n")


def get_nsw(file_path):
    with open(file_path, "r") as file:
        nsw_list = []
        for line in file:
            if re.search("Ionic steps remaining", line):
                nsw_list.append(int(line.split()[-1]))
    return nsw_list[-1]


if __name__ == '__main__':
    args = read_command_line_arguments()

    atoms_1 = None

    if args.restart:  # ase.io should read the structure of last index in vasprun_count.xml
        vasprun_files = [file for file in os.listdir() if "vasprun_" in file]
        if vasprun_files:
            vasprun_files.sort(key=lambda x: int(re.split("[_.]", x)[1]))
            args.file = vasprun_files[-1]
            args.index = -1
            args.nsw = get_nsw(file_path="./job.out") - \
                len(get_atoms(file=args.file, index=':'))

    atoms_1 = get_atoms(file=args.file, index=args.index)
    print_verbose(args)

    if args.new_lattice:
        a = args.new_lattice
        atoms_1.set_cell(cell=[a, a, a])

    run_isif_2(atoms_1, kwargs=args)
