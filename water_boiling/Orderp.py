# -*- coding: utf-8 -*-
# Copyright (c) 2015, PyRETIS Development Team.
# Distributed under the LGPLv3 License. See LICENSE for more info.
"""This file defines the order parameter used for the hydrate example."""
import logging
import os
import subprocess
import numpy as np
import mdtraj as md
from infretis.classes.orderparameter import OrderParameter
logger = logging.getLogger(__name__)  # pylint: disable=C0103
logger.addHandler(logging.NullHandler())


class H20Hole(OrderParameter):
    """H20Ccontinuous(OrderParameter).

    This class counts the size of H2O holes

    Attributes
    ----------
    name : string
        A human readable name for the order parameter
    """

    def __init__(self):
        """Initialize the order parameter.

        """
        super().__init__(description='H2O-hole')
        # trj=md.load('gromacs_input/conf.gro')  # <---- shortcut, make sure the file is there
        # self.idx_o = trj.top.select("symbol == O")


    def calculate(self, system):
        """Calculate the order parameter.

        Here, the order parameter is just the minimal value of 
        local H2O density in the thin film (2D density).

        Parameters
        ----------
        system : object like :py:class:`.System`
            This object is used for the actual calculation, typically
            only `system.particles.pos` and/or `system.particles.vel`
            will be used. In some cases `system.forcefield` can also be
            used to include specific energies for the order parameter.

        Returns
        -------
        out : float
            The order parameter list.
        """

        dens = (len(system.pos)/3)/system.box[0]**3

        return [-dens*(18.0153/((10**-21)*6.022*10**23))]
