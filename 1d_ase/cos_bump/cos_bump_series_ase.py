"""ASE-compatible cosine bump series potential.

This module defines a 1D potential with series of cosine bump barriers
and optional harmonic walls at the boundaries.
Converted from PyRETIS version, EW 2025
"""
import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.constraints import FixCartesian


class CosBumpSeriesWalls1D(Calculator):
    """ASE-compatible 1D cosine bump series potential with walls.
    
    The potential consists of:
    - Harmonic walls outside [xleft, xright]
    - Flat regions between walls and bumps
    - Series of cosine bumps between bleft and bright
    
    V(x) = k * (x - xleft)^2 / 2       if x < xleft
    V(x) = 0                            if xleft <= x < bleft
    V(x) = deltaf/2 * (1 - cos(...))   if bleft <= x <= bright2
    V(x) = 0                            if bright2 < x <= xright
    V(x) = k * (x - xright)^2 / 2      if x > xright
    
    where bright2 = bleft + nbump * (bright - bleft)
    
    Parameters:
        xleft (float): Left boundary for walls (default: -0.5)
        xright (float): Right boundary for walls (default: 0.5)
        k (float): Spring constant for the walls (default: 0.0)
        bleft (float): Left edge of first bump (default: -0.1)
        bright (float): Right edge of first bump (default: 0.1)
        deltaf (float): Height of bumps in units of kBT (default: 1.0)
        nbump (int): Number of bumps (default: 1)
        auto_constrain (bool): Auto-apply 1D constraint (default: True)
    """
    
    implemented_properties = ['energy', 'forces']

    def __init__(self, xleft=-0.5, xright=0.5, k=0.0, 
                 bleft=-0.1, bright=0.1, deltaf=1.0, nbump=1,
                 auto_constrain=True, **kwargs):
        super().__init__(**kwargs)
        self.xleft = xleft
        self.xright = xright
        self.k = k
        self.bleft = bleft
        self.bright = bright
        self.deltaf = deltaf
        self.nbump = nbump
        self.auto_constrain = bool(auto_constrain)
        self._constraint_applied = False
        
        # Validate parameters
        if self.bleft >= self.bright:
            import warnings
            warnings.warn("bleft should be < bright!")
        if self.bleft < self.xleft:
            import warnings
            warnings.warn("bleft should be >= xleft!")
        if self.bright > self.xright:
            import warnings
            warnings.warn("bright should be <= xright!")
        if self.xleft >= self.xright:
            import warnings
            warnings.warn("xleft should be < xright!")

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        
        # Auto-apply constraint (fix y/z) the first time calculate is called
        if getattr(self, 'auto_constrain', False) and not getattr(self, '_constraint_applied', False):
            atoms.set_constraint(create_1d_constraint(atoms))
            self._constraint_applied = True

        positions = atoms.get_positions()[:, 0]  # 1D x-coordinates
        
        # Derived quantities
        L = self.bright - self.bleft  # width of one bump
        bright2 = self.bleft + self.nbump * L  # end of all bumps
        xright2 = bright2 + (self.xright - self.bright)  # adjusted right wall
        
        energy = 0.0
        forces = np.zeros_like(atoms.get_positions())

        for i, x in enumerate(positions):
            # Energy calculation
            if x < self.xleft:
                # Left wall: V = k * (x - xleft)^2 / 2
                energy += 0.5 * self.k * (x - self.xleft)**2
                # Force: F = -dV/dx = -k * (x - xleft)
                forces[i, 0] = -self.k * (x - self.xleft)
                
            elif x > xright2:
                # Right wall: V = k * (x - xright2)^2 / 2
                energy += 0.5 * self.k * (x - xright2)**2
                # Force: F = -dV/dx = -k * (x - xright2)
                forces[i, 0] = -self.k * (x - xright2)
                
            elif (x < self.bleft) or (x > bright2):
                # Flat regions: V = 0, F = 0
                forces[i, 0] = 0.0
                
            else:
                # Cosine bump region
                # V = deltaf/2 * (1 - cos(2*pi*(x-bleft)/L))
                energy += self.deltaf / 2.0 * (1.0 - np.cos(2.0 * np.pi * (x - self.bleft) / L))
                # F = -dV/dx = -deltaf*pi/L * sin(2*pi*(x-bleft)/L)
                forces[i, 0] = -self.deltaf * np.pi / L * np.sin(2.0 * np.pi * (x - self.bleft) / L)

        self.results['energy'] = energy
        self.results['forces'] = forces

    def set_parameters(self, **params):
        """Update potential parameters.
        
        Parameters:
            xleft, xright (float): Wall boundaries
            k (float): Spring constant
            bleft, bright (float): Bump region boundaries
            deltaf (float): Bump height
            nbump (int): Number of bumps
        """
        for key, value in params.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                import warnings
                warnings.warn(f'Ignored unknown parameter "{key}"')
        
        # Validate updated parameters
        if self.bleft >= self.bright:
            import warnings
            warnings.warn("bleft should be < bright!")
        if self.xleft >= self.xright:
            import warnings
            warnings.warn("xleft should be < xright!")

    def get_boundaries(self):
        """Return the boundaries of the potential regions.
        
        Returns:
            dict: Dictionary with xleft, xright, bleft, bright
        """
        return {
            'xleft': self.xleft,
            'xright': self.xright,
            'bleft': self.bleft,
            'bright': self.bright
        }

    def get_parameters(self):
        """Return all potential parameters.
        
        Returns:
            dict: All parameters
        """
        return {
            'xleft': self.xleft,
            'xright': self.xright,
            'k': self.k,
            'bleft': self.bleft,
            'bright': self.bright,
            'deltaf': self.deltaf,
            'nbump': self.nbump
        }

    def potential_at_position(self, x):
        """Calculate potential energy at a single position.
        
        Parameters:
            x (float): Position coordinate
            
        Returns:
            float: Potential energy at position x
        """
        L = self.bright - self.bleft
        bright2 = self.bleft + self.nbump * L
        xright2 = bright2 + (self.xright - self.bright)
        
        if x < self.xleft:
            return 0.5 * self.k * (x - self.xleft)**2
        elif x > xright2:
            return 0.5 * self.k * (x - xright2)**2
        elif (x < self.bleft) or (x > bright2):
            return 0.0
        else:
            return self.deltaf / 2.0 * (1.0 - np.cos(2.0 * np.pi * (x - self.bleft) / L))

    def force_at_position(self, x):
        """Calculate force at a single position.
        
        Parameters:
            x (float): Position coordinate
            
        Returns:
            float: Force at position x
        """
        L = self.bright - self.bleft
        bright2 = self.bleft + self.nbump * L
        xright2 = bright2 + (self.xright - self.bright)
        
        if x < self.xleft:
            return -self.k * (x - self.xleft)
        elif x > xright2:
            return -self.k * (x - xright2)
        elif (x < self.bleft) or (x > bright2):
            return 0.0
        else:
            return -self.deltaf * np.pi / L * np.sin(2.0 * np.pi * (x - self.bleft) / L)


def create_1d_constraint(atoms):
    """Return a FixCartesian constraint that fixes motion to the x-axis.

    Usage:
        atoms.set_constraint(create_1d_constraint(atoms))
    """
    # Fix y and z (allow x motion only) for each atom
    return [FixCartesian(i, [False, True, True]) for i in range(len(atoms))]
