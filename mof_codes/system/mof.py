from typing import Union, Optional
import ase
import ase.io
from base_class import System


class MOF(System):
    """Class to extract MOF info, manipulate linkers, and metal nodes, & also can add adsorbates"""

    def __init__(self, atoms: Union[ase.atoms.Atoms, str], **kwargs: Union[int, list]) -> None:
        """Will create an instance/ object of class MOF.
        
        Parameters
        ----------
        atoms: ase.atoms.Atoms or str
            Can either be an ase.atoms.Atoms object or as a str containing the file path of the MOF structure.
        kwargs: int, list
            Can take any number of keyword arguments.
        """

        super().__init__(atoms, **kwargs)
        self.metal_atoms: Optional[list[ase.atoms.Atom]] = None
        self.linker_atoms: Optional[list[ase.atoms.Atom]] = None

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


if __name__ == '__main__':
    mof = MOF("./ase_atoms/other_structures/ZIF-8.traj", metal_elements=['Zn'])
    print(mof._get_index_connectivity(mof.atoms[0]))
    print(type(mof.linker_atoms[0]))
