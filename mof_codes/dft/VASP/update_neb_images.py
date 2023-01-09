#!/usr/bin/env python

import ase.io
from ase.parallel import *
import os
import argparse
import numpy as np


def write_pre_neb_images(dirs):
    all_atoms: list(list(ase.atoms.Atoms)) = []
    for dir in dirs[1:-1]:  # not reading IS and FS folders
        os.chdir(f"./{dir}")
        os.system(
            'mkdir -p temp; cp CONTCAR.pre-neb ./temp/CONTCAR; cp OUTCAR.pre-neb ./temp/OUTCAR; cp OSZICAR.pre-neb ./temp/OSZICAR')
        all_atoms.append(ase.io.read(f"./temp/OUTCAR", index=":"))
        os.system('rm -rf temp')
        os.chdir("..")

    all_atoms_flat: list(ase.atoms.Atoms) = [
        atoms for image in all_atoms for atoms in image]
    nsw_pre_neb = len(all_atoms_flat) // len(dirs[1:-1])

    with open("./neb_stats.txt", "w") as fi:
        i = 0
        traj = ase.io.Trajectory("./all_neb_images.traj", "w")
        fi.write(f"Pre-NEB details (lclimb is False)\n\n")
        while(i < nsw_pre_neb):
            fi.write(f"NEB {i+1}:\n")
            for j, dir in enumerate(dirs[1:-1]):
                atoms = all_atoms[j][i]
                f = atoms.get_forces()
                energy = atoms.get_potential_energy()
                forces = np.sqrt(
                    np.square(f[:, 0]) + np.square(f[:, 1]) + np.square(f[:, 2]))
                os.system(
                    f"grep -w \"{i+1} F\" {dir}/OSZICAR.pre-neb" + " | awk '{print $NF}' > temp")
                f = open("temp", "r")
                magmom = float(f.read())
                f.close()
                fi.write(
                    f"{dir} \t {energy:.4f} \t {max(forces):.4f} \t {magmom:.4f} \n")
                traj.write(atoms)
            fi.write("\n")
            i += 1
        traj.close()
        os.system('rm temp')


def write_actual_neb_images(dirs):
    all_atoms: list(list(ase.atoms.Atoms)) = []
    for dir in dirs[1:-1]:  # not reading IS and FS folders
        all_atoms.append(ase.io.read(f"./{dir}/OUTCAR", index=":"))
    min_nsw = min([len(atoms_list) for atoms_list in all_atoms])

    if min_nsw > 1:
        with open("./neb_stats.txt", "a") as fi:
            neb_count = None

            """below lines finds correct NEB iteration number immediately after pre-neb."""
            os.system('grep \"Actual-NEB\" neb_stats.txt > temp')
            if os.stat("./temp").st_size == 0:
                os.system(
                    "grep \"NEB\" neb_stats.txt | tail -n1 | awk '{print $NF}' > temp1")
                f = open("temp1", "r")
                neb_count = int(f.read().replace(":", ""))
                f.close()
                os.system('rm temp1')
                fi.write(f"Actual-NEB details (lclimb is True)\n\n")
            os.system('rm temp')

            """below lines finds correct NEB iteration number once actual neb has started."""
            if not neb_count:
                os.system(
                    "grep \"NEB\" neb_stats.txt | tail -n1 | awk '{print $NF}' > temp1")
                f = open("temp1", "r")
                neb_count = int(f.read().replace(":", ""))
                f.close()
                os.system('rm temp1')

            i = 0
            traj = ase.io.Trajectory("./all_neb_images.traj", "a")
            while(i < min_nsw-1):  # just to ensure ediff tolerance (increasing the chances)
                fi.write(f"NEB {neb_count + i + 1}:\n")
                for j, dir in enumerate(dirs[1:-1]):
                    atoms = all_atoms[j][i]
                    f = atoms.get_forces()
                    energy = atoms.get_potential_energy()
                    forces = np.sqrt(
                        np.square(f[:, 0]) + np.square(f[:, 1]) + np.square(f[:, 2]))
                    os.system(
                        f"grep -w \"{i+1} F\" {dir}/OSZICAR" + " | awk '{print $NF}' > temp")
                    f = open("temp", "r")
                    magmom = float(f.read())
                    f.close()
                    fi.write(
                        f"{dir} \t {energy:.4f} \t {max(forces):.4f} \t {magmom:.4f} \n")
                    traj.write(atoms)
                fi.write("\n")
                i += 1
            traj.close()
            os.system('rm temp')


def write_path_images():
    """This function will write/append the structures of the path images from their
    corresponding folder."""

    dirs = [name for name in os.listdir(
        ".") if os.path.isdir(name) and name[0] == "0"]

    if not dirs:
        raise "No folders found of NEB images with standard naming convention."

    dirs.sort(key=lambda x: int(x[1]))  # dirs = ["00", "01", "02", ...]

    if not os.path.exists("./all_neb_images.traj"):
        write_pre_neb_images(dirs)

    if os.path.exists("./all_neb_images.traj"):
        write_actual_neb_images(dirs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=f"Appends the neb images in \"all_neb_images.traj\" inside the NEB calculation directory")
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
        write_path_images()
        os.chdir(start_dir)
