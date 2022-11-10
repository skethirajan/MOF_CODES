from typing import Union, Optional
import ase
import ase.io
import ase.cell
from ase.neighborlist import natural_cutoffs, NeighborList
from scipy.sparse.csgraph import connected_components
import numpy as np


class System:
    """Base Class that children classes can inherit."""

    @staticmethod
    def _read_atoms(atoms: str, index: int) -> ase.atoms.Atoms:
        """Method to read the given structure file.

        Parameters
        ----------
        atoms: str
            Path of the file containing the given structure.
        index: int
            Structure at this configuration will be used.

        Returns
        -------
        ase.atoms.Atoms object of the given structure
        """

        if isinstance(atoms, str):
            if ".traj" in atoms:
                return ase.io.read(filename=atoms, index=index, format="traj")
            elif ".xyz" in atoms:
                return ase.io.read(filename=atoms, index=index, format="xyz")
            elif ".xml" in atoms:
                return ase.io.read(filename=atoms, index=index, format="vasp-xml")
            else:
                raise f"Check the format of the file, {atoms} to read."

    @staticmethod
    def _convert_to_ase_atoms(atoms: list[ase.atoms.Atom], cell: ase.cell.Cell) -> ase.atoms.Atoms:
        return ase.atoms.Atoms(atoms, cell=cell, pbc=True)

    @staticmethod
    def _mapping_indices(atoms: list[ase.atoms.Atom]) -> dict[int, int]:
        map_indices = {}  # mapping from index of ase.atoms.Atoms i.e. self.atoms to index of list[ase.atoms.Atom]
        for ind, atom in enumerate(atoms):
            map_indices.update({ind: atom.index})
        return map_indices

    def __init__(self, atoms: Union[ase.atoms.Atoms, list[ase.atoms.Atom], str],
                 **kwargs: Union[int, list, ase.cell.Cell]) -> None:
        """Will create an instance/ object of class.

        Parameters
        ----------
        atoms: ase.atoms.Atoms or str or list of ase.atoms.Atom
            Can either be an ase.atoms.Atoms object or as a str containing the file path of the given structure. It can
            also be as list[ase.atoms.Atom], say extracting linker atoms from a MOF.
        kwargs: int, list, ase.cell.Cell
            Can take any number of keyword arguments (optional).
        """

        self.atoms: Optional[ase.atoms.Atoms] = None
        self.is_ase_atoms: bool = True

        # becomes dict[int, int] only when self.is_ase_atoms is false
        self.map_indices: Optional[dict[int, int]] = None

        # becomes list[ase.atoms.Atom] only when self.is_ase_atoms is false
        self.list_of_atom: Optional[list[ase.atoms.Atom]] = None

        if isinstance(atoms, str):  # atoms is a file
            index = kwargs["index"] if kwargs.get("index") else -1
            self.atoms = self._read_atoms(atoms, index)
        elif isinstance(atoms, ase.atoms.Atoms):  # atoms is an ase.atoms.Atoms type
            self.atoms = atoms
        elif isinstance(atoms, list):  # atoms is a list type
            if all(isinstance(i, ase.atoms.Atom) for i in atoms):  # all elements in atoms list is ase.atoms.Atom type
                if kwargs.get("cell") and isinstance(kwargs["cell"], ase.cell.Cell):
                    self.is_ase_atoms = False
                    self.atoms = self._convert_to_ase_atoms(atoms, kwargs["cell"])
                    self.map_indices = self._mapping_indices(atoms)
                    self.list_of_atom = atoms

        if not self.atoms:
            raise "self.atoms is None. Check the atoms argument while creating an instance of the Class."

    @property
    def _get_natural_cutoffs(self) -> list[float]:
        """Method to calculate natural cutoffs of atoms in the given structure based on covalent radii.

        Returns
        -------
        ase atoms length list of float values based on the covalent radii.
        """
        return natural_cutoffs(self.atoms)

    @property
    def _get_neighborlist(self) -> ase.neighborlist.NeighborList:
        """Method that returns an instance of the ase.neighborlist.NeighborList class based on natural_cutoffs.

        Returns
        -------
        ase.neighborlist.NeighborList class object for "atoms" attribute of Class.
        """

        neighborlist = NeighborList(cutoffs=self._get_natural_cutoffs, self_interaction=False, bothways=True)
        neighborlist.update(self.atoms)
        return neighborlist

    @property
    def _connected_components(self) -> tuple[int, np.ndarray]:
        nl = self._get_neighborlist
        return connected_components(nl.get_connectivity_matrix())

    def _get_index_connectivity(self, atom: ase.atoms.Atom) -> list[int]:
        """Method that returns indices of first neighbor atoms as list for a given atom.

        Parameters
        ----------
        atom: ase.atoms.Atom
            Atom whose first neighbors to be found.

        Returns
        -------
        list of indices of primary neighbors for the given atom.
        """
        neighborlist = self._get_neighborlist
        neighbor_indices = list(neighborlist.get_neighbors(atom.index)[0])
        if self.is_ase_atoms:
            return neighbor_indices  # these indices are correct
        else:
            return [self.map_indices[ind] for ind in neighbor_indices]  # return indices w.r.t list[ase.atoms.Atom]

    @property
    def _get_indices_connectivity(self) -> dict[int, list[int]]:
        """Method that returns indices of first neighbor atoms for all the atoms as dictionary for the given structure.

        Returns
        -------
        Dictionary with key as atom index and values as a list of its neighbor indices.
        """
        indices_connectivity = {}
        for atom in self.atoms:
            neighbor_indices = self._get_index_connectivity(atom)
            if self.is_ase_atoms:
                indices_connectivity.update({atom.index: neighbor_indices})
            else:
                indices_connectivity.update({self.map_indices[atom.index]: neighbor_indices})
        return indices_connectivity


# if __name__ == '__main__':
#     system = System(atoms="./ase_atoms/other_structures/ZIF-8.traj", metal_elements=['Zn'])
#     print('0:', system._get_index_connectivity(system.atoms[0]))
#     print('1:', system._get_index_connectivity(system.atoms[1]))
#     # atoms1 = ase.io.read("./ase_atoms/other_structures/ZIF-8.traj", format="traj", index=-1)
#     # system = System(atoms=atoms1)
#     # print(system._get_index_connectivity(system.atoms[0]))
#     atoms_ = ase.io.read("./ase_atoms/other_structures/ZIF-8.traj", format="traj", index=-1)
#     atoms1 = [atoms_[0], atoms_[1], atoms_[216], atoms_[238], atoms_[2], atoms_[96], atoms_[46]]
#     cell1 = atoms_.cell
#     system = System(atoms=atoms1, cell=cell1)
#     # print(system.map_indices)
#     # print(system._get_index_connectivity(system.atoms[0]))
#     print(system._get_indices_connectivity)
