#!/usr/bin/env python

import ase.io
from ase.parallel import *
import os
import argparse
import numpy as np


def write_converged_images():
    """This function will append the converged path of the images from their
    corresponding folder."""

    os.system('grep \"accuracy\" vasp.out > temp')

    if os.stat("./temp").st_size == 0:
        raise "Not reached required accuracy i.e. NEB is not yet converged!!!"

    else:
        dirs = [name for name in os.listdir(
            ".") if os.path.isdir(name) and name[0] == "0"]

        if not dirs:
            raise "No folders found of NEB images with standard naming convention."

        dirs.sort(key=lambda x: int(x[1]))  # dirs = ["00", "01", "02", ...]

        all_atoms: list(ase.atoms.Atoms) = []
        for dir in dirs[1:-1]:  # not reading IS and FS folders
            all_atoms.append(ase.io.read(f"./{dir}/OUTCAR", index=-1))

        with open("./neb_stats.txt", "a") as fi:
            os.system(
                "grep \"NEB\" neb_stats.txt | tail -n1 | awk '{print $NF}' > temp1")
            f = open("temp1", "r")
            neb_count = int(f.read().replace(":", ""))
            f.close()
            os.system('rm temp1')
            traj = ase.io.Trajectory("./all_neb_images.traj", "a")

            fi.write(f"NEB {neb_count + 1}:\n")
            for j, dir in enumerate(dirs[1:-1]):
                atoms = all_atoms[j]
                f = atoms.get_forces()
                energy = atoms.get_potential_energy()
                forces = np.sqrt(
                    np.square(f[:, 0]) + np.square(f[:, 1]) + np.square(f[:, 2]))
                os.system(
                    f"grep \"F\" {dir}/OSZICAR | tail -n 1" + " | awk '{print $NF}' > temp")
                f = open("temp", "r")
                magmom = float(f.read())
                f.close()
                fi.write(
                    f"{dir} \t {energy:.4f} \t {max(forces):.4f} \t {magmom:.4f} \n")
                traj.write(atoms)
            traj.close()
            os.system('rm temp')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=f"Appends the converged neb images in \"all_neb_images.traj\" inside the NEB calculation directory")
    parser.add_argument("--calc_dir", type=str, default=".",
                        help="path of NEB calculation directory (default is \".\" i.e. current working directory)", metavar="")

    args = parser.parse_args()
    if args.calc_dir[-1] == "/":
        args.calc_dir = args.calc_dir[:-1]

    if not os.path.exists(args.calc_dir):
        raise f"Given path of calculation directory does not exist."
    else:
        start_dir = os.getcwd()
        os.chdir(args.calc_dir)
        write_converged_images()
        os.chdir(start_dir)
