import numpy as np
import ase.io
from ase.atoms import Atoms
from ase.neighborlist import mic
from ase.visualize import view


class AddAdsorbate:
    def __init__(self, system_atoms, adsorbate_atoms=Atoms(), adsorbates_to_add=1, max_iter=10000, min_distance=2.0,
                 set_seed=0, start_tag=None):
        self.system_initial = system_atoms
        self.cell_arr = self.system_initial.cell.array
        self.adsorbate_initial = adsorbate_atoms

        self.adsorbates_to_add = adsorbates_to_add
        self.adsorbates_added = 0
        self.max_iter = max_iter
        self.min_distance = min_distance
        self.tag = start_tag

        self.system = self.system_initial.copy()
        self.adsorbate = self.adsorbate_initial.copy()

        np.random.seed(set_seed)

    def rotate_adsorbate(self):
        u, v = np.random.normal(size=3), np.random.normal(size=3)
        u, v = u / np.linalg.norm(u), v / np.linalg.norm(v)
        self.adsorbate.rotate(u, v)

    def move_adsorbate(self):
        cop = np.sum(self.adsorbate.positions, axis=0) / len(self.adsorbate)
        distance_to_move = self.random_position - cop
        self.adsorbate.translate(displacement=distance_to_move)

    @property
    def random_position(self):
        random_positions_inside_cell = self.cell_arr * \
            np.random.random(size=(3, 1))  # this is a 3 X 3 array
        # this is a point inside the cell
        random_position = np.sum(random_positions_inside_cell, axis=0)
        return random_position

    @property
    def adsorbate_distance(self):
        distances_matrix = np.zeros(
            shape=(len(self.adsorbate), len(self.system)))
        for i, ads_atom in enumerate(self.adsorbate):
            distance_arr = ads_atom.position - \
                self.system.positions  # this is len(system) X 3
            distance = mic(dr=distance_arr, cell=self.cell_arr,
                           pbc=True)  # this is len(system) X 3
            distances_matrix[i] = np.linalg.norm(distance, axis=1)
        return distances_matrix

    def add_adsorbate(self):
        self.rotate_adsorbate()
        self.move_adsorbate()

        iter_count = 0

        while np.min(self.adsorbate_distance) <= self.min_distance:
            if iter_count == self.max_iter:
                return False
            else:
                self.rotate_adsorbate()
                self.move_adsorbate()
                iter_count += 1
        return True

    def fill_system(self, **kwargs):
        if kwargs.get("adsorbate_atoms"):
            self.adsorbate_initial = kwargs["adsorbate_atoms"]

        if kwargs.get("adsorbates_to_add"):
            self.adsorbates_to_add = kwargs["adsorbates_to_add"]

        if kwargs.get("max_iter"):
            self.max_iter = kwargs["max_iter"]

        if kwargs.get("min_distance"):
            self.min_distance = kwargs["min_distance"]

        if kwargs.get("start_tag"):
            self.tag = kwargs["start_tag"]

        while self.adsorbates_added < self.adsorbates_to_add:
            self.adsorbate = self.adsorbate_initial.copy()
            if self.add_adsorbate():
                if self.tag:
                    self.adsorbate.set_tags(self.tag)
                    self.tag += 1
                self.system += self.adsorbate
                self.adsorbates_added += 1
            else:
                raise f"Unable to add adsorbate {self.adsorbates_added + 1} to the system."


"""
if __name__ == "__main__":
    linkers = ['linker-1(1)_relax.traj', 'linker-2(1)_relax.traj', 'linker-3(1)_relax.traj',
               'linker-4(1)_relax.traj', 'linker-5(1)_relax.traj', 'linker-6(1)_relax.traj']
    mofs = ['mof-1_relax.traj', 'mof-2_relax.traj', 'mof-3_relax.traj',
            'mof-4_relax.traj', 'mof-5_relax.traj', 'mof-6_relax.traj']

    for file_1, file_2 in zip(mofs, linkers):
        mof_file = "./mofs_1/relaxed_structures/mofs/" + file_1
        ads_file = "./mofs_1/relaxed_structures/linkers/" + file_2

        adsorbate = ase.io.read(ads_file, format="traj", index=-1)
        mof = ase.io.read(mof_file, format="traj", index=-1)

        mof_name = file_1.split(sep="_")[0]
        adsorbate_name = file_2.split(sep="_")[0].replace("(1)", "")

        for count in [1, 2, 4]:
            ads = AddAdsorbate(
                system_atoms=mof, adsorbate_atoms=adsorbate, adsorbates_to_add=count, start_tag=1)
            ads.fill_system()
            file_name = f"./mofs_1/initial_structures/{mof_name}_{adsorbate_name}({str(count)})_initial.traj"
            ase.io.write(file_name, ads.system)
            print(f"{file_name} created!")
"""
