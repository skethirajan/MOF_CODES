import argparse
import ase.io
from ase.neb import NEB
import os


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
    parser = argparse.ArgumentParser(
        description="Interpolates between IS to TS and TS to FS")
    parser.add_argument("--IS", type=str, required=True,
                        help="path of initial state structure", metavar="")
    parser.add_argument("--TS", type=str, required=True,
                        help="path of guess transition state structure", metavar="")
    parser.add_argument("--FS", type=str, required=True,
                        help="path of final state structure", metavar="")
    parser.add_argument("--out_dir", type=str, default=".",
                        help="output directory (default is \".\" i.e. current working directory)", metavar="")
    parser.add_argument("--IS_index", type=int, default=-1,
                        help="index of IS structure (default is last index)", metavar="")
    parser.add_argument("--TS_index", type=int, default=-1,
                        help="index of TS structure (default is last index)", metavar="")
    parser.add_argument("--FS_index", type=int, default=-1,
                        help="index of FS structure (default is last index)", metavar="")
    parser.add_argument("--method", type=str, default="linear", choices=[
                        "linear", "idpp"], help="NEB interpolation method (default is linear)", metavar="")
    parser.add_argument("--left", type=int, default=2,
                        help="number of images to interpolate between IS and TS (default is 2)", metavar="")
    parser.add_argument("--right", type=int, default=2,
                        help="number of images to interpolate between TS and FS (default is 2)", metavar="")
    args = parser.parse_args()

    if args.out_dir[-1] == "/":
        args.out_dir = args.out_dir[:-1]

    return args


if __name__ == '__main__':
    args = read_command_line_arguments()

    initial_structure = get_atoms(file=args.IS, index=args.IS_index)
    initial_structure.pbc = True
    transition_structure = get_atoms(file=args.TS, index=args.TS_index)
    transition_structure.pbc = True
    final_structure = get_atoms(file=args.FS, index=args.FS_index)
    final_structure.pbc = True

    guess_path = []
    images = [initial_structure]
    images += [initial_structure.copy() for _ in range(args.left)]
    images += [transition_structure]
    neb = NEB(images)
    neb.interpolate(method=args.method, mic=True)
    guess_path = images[1:]  # this includes the guess TS but excludes IS
    del images

    images = [transition_structure]
    images += [transition_structure.copy() for i in range(args.right)]
    images += [final_structure]
    neb = NEB(images)
    neb.interpolate(method=args.method, mic=True)
    guess_path += images[1:-1]  # this excludes both the guess TS and FS
    del images

    if not os.path.exists(f"{args.out_dir}") and args.out_dir != ".":
        os.makedirs(f"{args.out_dir}")
    ase.io.write(f"{args.out_dir}/is.traj", initial_structure)
    ase.io.write(f"{args.out_dir}/fs.traj", final_structure)
    ase.io.write(f"{args.out_dir}/guess_path.traj", guess_path)
