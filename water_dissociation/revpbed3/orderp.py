import numpy as np
import logging
from MDAnalysis.analysis import distances
from infretis.classes.orderparameter import OrderParameter
import ase
from ase.io import read

logger = logging.getLogger(__name__)  # pylint: disable=invalid-name
logger.addHandler(logging.NullHandler())

# Assign global variables:
# BOX = [9.8528, 9.8528, 9.8528, 90, 90, 90]
# NATOMS = 32*3
# HIDX = [i for i in range(NATOMS) if i % 3 != 0]
# OIDX = list(range(0, 96, 3))


class OrderX(OrderParameter):
    """A positional order parameter.

    Order parameter for the hysteresis example. In addition to using
    the position, we also use the energy to tell if we are in states A/B.


    Attributes
    ----------
    index : integer
        This is the index of the atom which will be used, i.e.
        system.particles.pos[index] will be used.
    inter_a : float
        An interface such that we are in state A for postions < inter_a.
    inter_b : float
        An interface such that we are in state B for postions > inter_b.
    energy_a : float
        An energy such that we are in state A for potential energy < energy_a.
    energy_b : float
        An energy such that we are in state A for potential energy < energy_b.
    dim : integer
        This is the dimension of the coordinate to use.
        0, 1 or 2 for 'x', 'y' or 'z'.
    periodic : boolean
        This determines if periodic boundaries should be applied to
        the position or not.

    """

    def __init__(self, index, base, periodic=False):
        """Initialise the order parameter.

        Parameters
        ----------
        index : tuple of ints
            This is the indices of the atoms we will use the position of.
        periodic : boolean, optional
            This determines if periodic boundary conditions should be
            applied to the position.

        """
        pbc = 'Periodic' if periodic else 'Non-periodic'
        txt = '{} distance, particles {} and {}'.format(
            pbc,
            index[0],
            index[1]
        )
        super().__init__(description=txt, velocity=False)
        self.periodic = periodic
        self.index = index
        self.atoms = read(base)
        self.box = [self.atoms.cell[i][i] for i in [0,1,2]] + [90]*3
        self.hindex = [i for i, j in enumerate(self.atoms.numbers) if j == 1]
        self.oindex = [i for i, j in enumerate(self.atoms.numbers) if j == 8]


    def calculate(self, system):
        """Calculate the order parameter.

        Here, the order parameter is just the distance between two
        particles.

        Parameters
        ----------
        system : object like :py:class:`.System`
            The object containing the positions and box used for The
            calculation.

        Returns
        -------
        out : list of floats
            The rate-of-change of the distance order parameter.
        """
        print('cheese', self.atoms.pbc)
        # set up system variables
        pos = system.pos

        # key: OIDX, value: dict of HIDXs and dists
        dinfo = {i: {'hidxs': [], 'dists': []} for i in self.oindex}

        # calc distance pairs
        for idx in self.hindex:
            dists = distances.distance_array(pos[idx],
                                             pos[self.oindex],
                                             box=self.box)[0]
            o_close = self.oindex[np.argmin(dists)]
            oh_dist = dists[np.argmin(dists)]
            dinfo[o_close]['hidxs'].append(idx)
            dinfo[o_close]['dists'].append(oh_dist)

        # check do not have isolated O:
        assert len(dinfo) == len(self.oindex)

        # calculate num h per o list
        olist = [len(i['hidxs']) for i in dinfo.values()]

        # if no dissociation:
        if all(i == 2 for i in olist):
            orderp = max([max(i['dists']) for i in dinfo.values()])
            hattach = -1
            oattach_min = -1
            oattach_max = -1
            typep = 0

        # we have 30 2's, h3o+, oh-
        elif sum([i == 2 for i in olist]) == 30 and set(olist) == set((1,2,3)):
        # elif sum([i == 2 for i in olist]) == len(olist) - 2:
            o_min = self.oindex[np.argmin(olist)]
            o_max = self.oindex[np.argmax(olist)]
            dists = distances.distance_array(pos[o_min],
                                             pos[dinfo[o_max]['hidxs']],
                                             box=BOX)[0]
            orderp = min(dists)
            hattach = dinfo[o_max]['hidxs'][np.argmin(dists)]
            oattach_min = o_min
            oattach_max = o_max
            typep = 1
        else:
            o_mins = np.array(self.oindex)[np.array(olist)==1]
            o_maxs = np.array(self.oindex)[np.array(olist)==3]
            dists = []
            # for o_min, o_max in zip(o_mins, o_maxs):
            hattach = -1
            for o_min in o_mins:
                for o_max in o_maxs:
                    dists_i = distances.distance_array(pos[o_min],
                                                       pos[dinfo[o_max]['hidxs']],
                                                       box=BOX)[0]
                    if dists and min(dists_i) < min(dists):
                        hattach = dinfo[o_max]['hidxs'][np.argmin(dists_i)]
                    dists += list(dists_i)
                    oattach_min = o_min
                    oattach_max = o_max
            orderp = min(dists)
            typep = 2

        return [orderp, typep, oattach_min, oattach_max, hattach]

