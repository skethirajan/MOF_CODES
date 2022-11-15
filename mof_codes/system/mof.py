from typing import Union, Optional
import ase
import ase.io
import ase.geometry
from base_class import System
from ase.visualize import view
import numpy as np


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
        self.mof_atoms: Optional[list[ase.atoms.Atom]] = None
        self.adsorbates: Optional[list[ase.atoms.Atom]] = None
        self.metal_atoms: Optional[list[ase.atoms.Atom]] = None
        self.linker_atoms: Optional[list[ase.atoms.Atom]] = None
        self.linker_clusters: Optional[dict[int, list[ase.atoms.Atom]]] = None

        self.adsorbates, self.mof_atoms = self._filter_structure()

        if kwargs.get("metal_elements"):
            self.set_metal_atoms(metal_elements=kwargs["metal_elements"])
        elif kwargs.get("metal_indices"):
            self.set_metal_atoms(metal_indices=kwargs["metal_indices"])

        if self.metal_atoms:  # only when self.metal_atoms is not None, self.set_linker_atoms() method is called
            self.set_linker_atoms()
            if self.linker_atoms:
                self.linker_clusters: dict[int, list[ase.atoms.Atom]] = self._get_linker_clusters()

    def _filter_structure(self):
        clusters = self._clusters
        adsorbates = {}
        mof_atoms = []
        if len(clusters.keys()) > 1:
            i = 1
            for cls_atoms in list(clusters.values()):
                if len(cls_atoms) < 50:  # number of atoms in each adsorbate should be less than 50 else it is a MOF
                    adsorbates.update({i: cls_atoms})
                    i += 1
                else:
                    mof_atoms.append(cls_atoms)
        else:
            mof_atoms = list(clusters.values())
        return adsorbates, *mof_atoms

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
            metal_atoms = [atom for atom in self.mof_atoms if atom.symbol in metal_elements]
        elif metal_indices:
            metal_atoms = [atom for atom in self.mof_atoms if atom.index in metal_indices]

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
            self.linker_atoms: list[ase.atoms.Atom] = [atom for atom in self.mof_atoms
                                                       if atom.index not in metal_indices]
        else:
            raise "First set metal nodes."

    def _is_metal_atom(self, atom_1: ase.atoms.Atom) -> bool:
        for atom_2 in self.metal_atoms:
            if self._compare_two_atoms(atom1=atom_1, atom2=atom_2):
                return True
        else:
            return False

    def _is_linker_atom(self, atom_1: ase.atoms.Atom) -> bool:
        for cluster_id, cluster_atoms in self.linker_clusters.items():
            for atom_2 in cluster_atoms:
                if self._compare_two_atoms(atom1=atom_1, atom2=atom_2):
                    return True
        else:
            return False

    def _is_adsorbate_atom(self, atom_1: ase.atoms.Atom) -> bool:
        for ads_id, ads_atoms in self.adsorbates.items():
            for atom_2 in ads_atoms:
                if self._compare_two_atoms(atom1=atom_1, atom2=atom_2):
                    return True
        else:
            return False

    def _get_linker_clusters(self) -> dict[int, list[ase.atoms.Atom]]:
        return System(atoms=self.linker_atoms, cell=self.atoms.cell)._clusters

    def _get_cluster_id_of_linker_atom(self, linker_atom_1: ase.atoms.Atom) -> Optional[int]:
        for cluster_id, cluster_atoms in self.linker_clusters.items():
            for atom_2 in cluster_atoms:
                if self._compare_two_atoms(atom1=linker_atom_1, atom2=atom_2):
                    return cluster_id
        else:
            return None

    def _get_cluster_ids_of_metal_atom(self, metal_atom_1: ase.atoms.Atom) -> Optional[list[int]]:
        neighbor_indices = self._get_index_connectivity(atom=metal_atom_1)
        neighbor_atoms = [self.atoms[index] for index in neighbor_indices]
        cluster_ids = []
        for atom_2 in neighbor_atoms:
            cluster_id = self._get_cluster_id_of_linker_atom(linker_atom_1=atom_2)
            if cluster_id:  # not None
                cluster_ids.append(cluster_id)
        return cluster_ids

    def _get_cluster_connectivity(self, cls_id: int, only_linkers: bool = True) \
            -> dict[ase.atoms.Atom, list[ase.atoms.Atom]]:
        cluster_atoms = self.linker_clusters[cls_id]
        metal_indices = [atom.index for atom in self.metal_atoms]
        cluster_connectivity = {}
        metal_indices_found = set()

        for atom1 in cluster_atoms:
            neighbor_indices = [index for index in self._get_index_connectivity(atom=atom1)]
            if only_linkers:
                allowed_indices = list(set(neighbor_indices) - set(metal_indices))
            else:
                allowed_indices = neighbor_indices
                if set(metal_indices) & set(allowed_indices):
                    metal_indices_found.update(set(metal_indices) & set(allowed_indices))
            allowed_indices.sort(key=lambda ind: self.atoms[ind].symbol)
            cluster_connectivity.update({atom1: [self.atoms[index] for index in allowed_indices]})

        if metal_indices_found:  # i.e. when only_linkers is false
            for index in metal_indices_found:
                metal_atom = self.atoms[index]
                allowed_metal_neighbors = []
                metal_neighbors = self._get_index_connectivity(atom=metal_atom)
                for i in metal_neighbors:
                    neighbor_atom = self.atoms[i]
                    if self._get_cluster_id_of_linker_atom(linker_atom_1=neighbor_atom) == cls_id:
                        allowed_metal_neighbors.append(neighbor_atom)
                else:
                    cluster_connectivity.update({metal_atom: allowed_metal_neighbors})
        return cluster_connectivity

    def extract_linker_molecule(self, cls_id, only_linkers=False) -> ase.atoms.Atoms:
        cluster_connectivity = self._get_cluster_connectivity(cls_id, only_linkers)
        atoms = ase.atoms.Atoms(list(cluster_connectivity.keys()), pbc=True, cell=self.atoms.cell)
        return self.center_molecule(atoms)

    @staticmethod
    def center_molecule(atoms: ase.atoms.Atoms) -> ase.atoms.Atoms:
        distances_with_mic = atoms.get_distances(a=0, indices=np.arange(0, len(atoms)), mic=True, vector=True)
        distances_without_mic = atoms.get_distances(a=0, indices=np.arange(0, len(atoms)), mic=False, vector=True)
        displacement = distances_with_mic - distances_without_mic
        atoms.translate(displacement)
        atoms.center()
        return atoms


if __name__ == '__main__':
    system = MOF("./ase_atoms/relaxed_structures/mof-6_relax.traj", metal_elements=['Zn'])
    molecule = system.extract_linker_molecule(cls_id=1, only_linkers=False)
    N_indices = [atom.index for atom in molecule if atom.symbol == 'N']
    Zn_indices = [atom.index for atom in molecule if atom.symbol == 'Zn']
    new_molecule = []
    for Zn_ind in Zn_indices:
        distances = list(molecule.get_distances(Zn_ind, indices=N_indices))
        N_index = distances.index(min(distances))
        molecule[Zn_ind].symbol = 'H'
        molecule.set_distance(a0=N_index, a1=Zn_ind, fix=0.5, distance=-1.6, indices=[Zn_ind], add=True)
    # view(molecule)
    # ase.io.write("./ase_atoms/initial_structures/linker-6_initial.traj", molecule)
    print(molecule)
