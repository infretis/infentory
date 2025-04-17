from infretis.classes.orderparameter import OrderParameter
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np

class PositionX(OrderParameter):
    """Position order parameter.


    Attributes
    ----------
    index : tuple of integers
        These are the indices used for the two particles.
        `system.particles.pos[index[0]]` and
        `system.particles.pos[index[1]]` will be used.
    periodic : boolean
        This determines if periodic boundaries should be applied to
        the distance or not.

    """

    def __init__(self, config):
        """Initialise order parameter.

        Parameters
        ----------
        index : tuple of ints
            This is the indices of the atom we will use the position of.
        periodic : boolean, optional
            This determines if periodic boundary conditions should be
            applied to the position.

        """
        super().__init__(description="dist", velocity=False)
        self.u = mda.Universe(config)
        self.pot = self.u.select_atoms("resname POT")
        self.cla = self.u.select_atoms("resname CLA")

    def calculate(self, system) -> list[float]:
        """Calculate the order parameter.

        Here, the order parameter is just the position
        of a particle.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The object containing the positions and box used for the
            calculation.

        Returns
        -------
        out : list of floats
            The distance order parameter.

        """
        # units ...
        self.u.atoms.positions = system.pos*10
        self.u.dimensions[:3] = system.box[:3]*10

        dist = distance_array(self.pot, self.cla, box=self.u.dimensions)[0][0]

        return [dist,]
