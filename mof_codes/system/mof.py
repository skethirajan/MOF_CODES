from typing import Union, Optional
import ase
import ase.io
from ase.neighborlist import natural_cutoffs, NeighborList


class MOF:
    """Class to extract MOF info, manipulate linkers, and metal nodes, & also can add adsorbates"""

    @staticmethod
    def read_atoms(atoms: str, index: int) -> ase.atoms.Atoms:
        """Method to read the MOF structure file.
        
        Parameters
        ----------
        atoms: str
            Path of the file containing MOF structure.
        index: int
            MOF structure at this configuration will be used.

        Returns
        -------
        ase.atoms.Atoms object of the MOF structure
        """

        if isinstance(atoms, str):
            if ".traj" in atoms:
                return ase.io.read(filename=atoms, index=index, format="traj")
            elif ".xyz" in atoms:
                return ase.io.read(filename=atoms, index=index, format="xyz")
            elif ".xml" in atoms:
                return ase.io.read(filename=atoms, index=index, format="vasp-xml")
            else:
                raise "Check the format of the file name to read."

    def __init__(self, atoms: Union[ase.atoms.Atoms, str], **kwargs: Union[int, list]) -> None:
        """Will create an instance/ object of class MOF.
        
        Parameters
        ----------
        atoms: ase.atoms.Atoms or str
            Can either be an ase.atoms.Atoms object or as a str containing the file path of the MOF structure.
        kwargs: int, list
            Can take any number of keyword arguments.
        """

        self.atoms: Optional[ase.atoms.Atoms] = None
        self.metal_atoms: Optional[list[ase.atoms.Atom]] = None
        self.linker_atoms: Optional[list[ase.atoms.Atom]] = None

        if not isinstance(atoms, ase.atoms.Atoms):
            index = kwargs["index"] if kwargs.get("index") else -1
            self.atoms = self.read_atoms(atoms, index)
        else:
            self.atoms = atoms

        if kwargs.get("metal_elements"):
            self.set_metal_atoms(metal_elements=kwargs["metal_elements"])
        elif kwargs.get("metal_indices"):
            self.set_metal_atoms(metal_indices=kwargs["metal_indices"])

        if self.metal_atoms:  # only when self.metal_atoms is not None, self.set_linker_atoms() method is called
            self.set_linker_atoms()

    def set_metal_atoms(self, metal_elements: Optional[list[str]] = None, metal_indices: Optional[list[int]] = None) \
            -> None:
        """Method to initialize the "metal_atoms" attribute of MOF class. One of the arguments should be assigned.

        Parameters
        ----------
        metal_elements: List of str
            Contains the elements of the metal nodes. Use this argument only when the elements of metal nodes are
            different from those in linkers and adsorbates if present.
        metal_indices: List of int
            Contains the ase.atoms.Atoms indices of the metal nodes. Use this argument when elements cannot be used.
        """

        metal_atoms = None
        if metal_elements:
            metal_atoms = [atom for atom in self.atoms if atom.symbol in metal_elements]
        elif metal_indices:
            metal_atoms = [atom for atom in self.atoms if atom.index in metal_indices]

        if metal_atoms:
            self.metal_atoms: list[ase.atoms.Atom] = metal_atoms
        else:
            raise "No metal nodes found with the given arguments."

    def set_linker_atoms(self) -> None:
        """Method to initialize the "linker_atoms" attribute of MOF class. "metal_atoms" attribute of MOF class should
        already be initialized to a non-Null value.
        """

        if self.metal_atoms:
            metal_indices = [atom.index for atom in self.metal_atoms]
            self.linker_atoms: list[ase.atoms.Atom] = [atom for atom in self.atoms if atom.index not in metal_indices]
        else:
            raise "First set metal nodes."

    @property
    def _get_natural_cutoffs(self) -> list[float]:
        """Method to calculate natural cutoffs of atoms in MOF based on covalent radii.

        Returns
        -------
        ase atoms length list of float values based on the covalent radii.
        """
        return natural_cutoffs(self.atoms)

    @property
    def _get_neighborlist(self) -> ase.neighborlist.NeighborList:
        """Method that returns an instance of the ase.neighborlist.NeighborList class based on natural_cutoffs of MOF.

        Returns
        -------
        ase.neighborlist.NeighborList class object for "atoms" attribute of MOF class.
        """

        neighborlist = NeighborList(cutoffs=self._get_natural_cutoffs, self_interaction=False, bothways=True)
        neighborlist.update(self.atoms)
        return neighborlist

    def _get_index_connectivity(self, atom: ase.atoms.Atom) -> list[int]:
        """Method that returns indices of first neighbor atoms as list for a given atom in MOF.

        Parameters
        ----------
        atom: ase.atoms.Atom
            Atom whose first neighbors to be found.

        Returns
        -------
        list of indices of primary neighbors for the given atom.
        """
        neighborlist = self._get_neighborlist
        return list(neighborlist.get_neighbors(atom.index)[0])

    @property
    def _get_indices_connectivity(self) -> dict[int, list[int]]:
        """Method that returns indices of first neighbor atoms for all the atoms as dictionary in the MOF.

        Returns
        -------
        Dictionary with key as atom index and values as a list of its neighbor indices.
        """
        indices_connectivity = {}
        for atom in self.atoms:
            neighbor_indices = self._get_index_connectivity(atom)
            indices_connectivity.update({atom.index: neighbor_indices})
        return indices_connectivity


if __name__ == '__main__':
    mof = MOF("./ase_atoms/other_structures/ZIF-8.traj", metal_elements=['Zn'])
    print(mof._get_index_connectivity(mof.atoms[0]))

