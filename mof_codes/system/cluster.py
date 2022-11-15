import ase.cell
from base_class import System


class Cluster(System):
    """Class to extract and manipulate linkers in a given MOF structure."""

    def __init__(self, atoms: list[ase.atoms.Atom], **kwargs: ase.cell.Cell):
        super().__init__(atoms, **kwargs)
        pass


if __name__ == '__main__':
    pass
