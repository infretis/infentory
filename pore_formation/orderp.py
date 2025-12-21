from infretis.classes.orderparameter import OrderParameter
import MDAnalysis as mda
from dztools.misc.mem_help import calc_chain
import numpy as np
from MDAnalysis.analysis.distances import distance_array


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
        super().__init__(description="chain", velocity=False)
        self.u = mda.Universe(config)
        self.config = config
        self.hoxy = "OH2 O11 O12 O13 O14"
        self.lip = "resname DMPC"
        self.coord_n = 26
        self.coord_d = 1.
        self.coord_r = 8.
        self.coord_z = 0.75
        self.lipid_p = self.u.select_atoms(f"name P")
        self.lmt_n = 12

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
        self.u.atoms.positions = system.pos * 10
        self.u.dimensions = list(system.box[:3]*10) + [90.]*3
        box = system.box


        # xi_ch, xi_p
        eps, eps_e, xcyl, ycyl, _ = calc_chain(self.u,
                                               lip = self.lip,
                                               hoxy = self.hoxy,
                                               coord_n = self.coord_n,
                                               coord_r = self.coord_r,
                                               coord_d = self.coord_d,
                                               coord_z = self.coord_z,
                                               )

        # d_llp
        xyz = self.lipid_p.atoms.positions
        radial = distance_array(np.array([xcyl, ycyl, 0]),
                                xyz*np.array([1, 1, 0]),
                                box = self.u.dimensions)[0]
        zavg = np.average(xyz[:, 2])
        ups, dws = [], []           # group idxes
        for i in np.argsort(radial):
            at = self.lipid_p.atoms[i]
            if xyz[i, 2] > zavg:
                ups.append(i)
            else:
                dws.append(i)
            if len(ups) >= 1 and len(dws) >= 1 and len(ups) + len(dws) >= self.lmt_n:
                break

        # get their positions, and get their pair stuff
        ups_pos = xyz[ups]
        dws_pos = xyz[dws]

        # d_llp
        d_llp = np.average(ups_pos[:, 2]) - np.average(dws_pos[:, 2])

        # linear regression line across xi_p, llp
        a = 0.25

        # op
        op = -d_llp/10 + (1/a) * eps_e

        # save op, xi_p, llp, xcyl, ycyl, box[0]
        return [op, eps_e, d_llp, xcyl, ycyl, box[0]]
