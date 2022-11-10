import ase
import ase.cell
from ase.neighborlist import natural_cutoffs, NeighborList
import numpy as np
from scipy.sparse.csgraph import connected_components


class Linker:
    """Class to extract and manipulate linkers in a given MOF structure."""

    def __init__(self, linker_atoms: list[ase.atoms.Atom], cell: ase.cell.Cell):
        self.atoms: list[ase.atoms.Atom] = linker_atoms
        self.cell: ase.cell.Cell = cell
        self.ase_atoms: ase.atoms.Atoms = self._convert_to_ase_atoms
        self.map_indices: dict[int, int] = self._mapping_indices
        self.clusters: dict[int, list[ase.atoms.Atom]] = self.get_clusters

    @property
    def _convert_to_ase_atoms(self) -> ase.atoms.Atoms:
        return ase.atoms.Atoms(self.atoms, cell=self.cell)

    @property
    def _mapping_indices(self) -> dict[int, int]:
        map_indices = {}  # mapping from index of ase.atoms.Atoms to index of list[ase.atoms.Atom]
        for ind, atom in enumerate(self.atoms):
            map_indices.update({ind: atom.index})
        return map_indices

    @property
    def _connected_components(self) -> tuple[int, np.ndarray]:
        nl = NeighborList(cutoffs=natural_cutoffs(self.ase_atoms), self_interaction=False, bothways=True)
        nl.update(self.ase_atoms)
        return connected_components(nl.get_connectivity_matrix())

    @property
    def get_clusters(self) -> dict[int, list[ase.atoms.Atom]]:
        linker_clusters = {}
        n_components, component_list = self._connected_components
        for cluster_id in range(0, n_components):
            indices = np.where(component_list == cluster_id)[0]  # gives indices according to ase.atoms.Atoms
            correct_indices = [self.map_indices[ind] for ind in indices]  # indices according to MOF.linker_atoms
            cluster = [self.atoms[index] for index in correct_indices]
            linker_clusters.update({cluster_id: cluster})
        return linker_clusters
