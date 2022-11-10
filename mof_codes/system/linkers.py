import ase
import ase.cell
import numpy as np
from mof import MOF
from base_class import System


class Linker(System):
    """Class to extract and manipulate linkers in a given MOF structure."""

    def __init__(self, atoms: list[ase.atoms.Atom], **kwargs: ase.cell.Cell):
        super().__init__(atoms, **kwargs)
        self.clusters: dict[int, list[ase.atoms.Atom]] = self._get_clusters

    @property
    def _get_clusters(self) -> dict[int, list[ase.atoms.Atom]]:
        linker_clusters = {}
        n_components, component_list = self._connected_components
        for cluster_id in range(0, n_components):
            indices = np.where(component_list == cluster_id)[0]  # gives indices according to ase.atoms.Atoms
            correct_indices = [self.map_indices[ind] for ind in indices]  # indices according to MOF.linker_atoms
            cluster = [self.atoms[index] for index in correct_indices]
            linker_clusters.update({cluster_id: cluster})
        return linker_clusters


if __name__ == '__main__':
    mof = MOF("./ase_atoms/other_structures/ZIF-8.traj", metal_elements=['Zn'])
    linkers = Linker(atoms=mof.linker_atoms, cell=mof.atoms.cell)
    print(linkers.clusters.keys())
