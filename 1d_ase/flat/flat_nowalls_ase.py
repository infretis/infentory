"""ASE-compatible flat walls potential.

This module defines a 1D flat walls potential for ASE simulations.
The potential is zero in the allowed region and has harmonic walls
at the boundaries.
EW, June 2025
"""
import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.constraints import FixCartesian


class FlatWalls1D(Calculator):
    """ASE-compatible 1D flat walls potential.
    
    The potential energy is zero within the allowed region [xleft, xright]
    and has harmonic walls outside:
    
    V(x) = 0                      if xleft <= x <= xright
    V(x) = k * (x - xleft)^2 / 2  if x < xleft  
    V(x) = k * (x - xright)^2 / 2 if x > xright
    
    Parameters:
        xleft (float): Left boundary of the flat region (default: -0.2)
        xright (float): Right boundary of the flat region (default: 1.0)
        k (float): Spring constant for the walls (default: 100.0)
    """
    
    implemented_properties = ['energy', 'forces']

    def __init__(self, xleft=-0.2, xright=0.4, k=0.0, auto_constrain=True, **kwargs):
        super().__init__(**kwargs)
        self.xleft = xleft
        self.xright = xright
        self.k = k
        self.auto_constrain = bool(auto_constrain)
        self._constraint_applied = False

        if self.xleft >= self.xright:
            import warnings
            warnings.warn("Setting xleft >= xright in FlatWalls1D potential!")

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        
        # Auto-apply constraint (fix y/z) the first time calculate is called
        if getattr(self, 'auto_constrain', False) and not getattr(self, '_constraint_applied', False):
            atoms.set_constraint(create_1d_constraint(atoms))
            self._constraint_applied = True

        positions = atoms.get_positions()[:, 0]  # 1D x-coordinates

        energy = 0.0
        forces = np.zeros_like(atoms.get_positions())

        for i, x in enumerate(positions):
            # Energy calculation
            if x < self.xleft:
                # Left wall: V = k * (x - xleft)^2 / 2
                energy += 0.5 * self.k * (x - self.xleft)**2
                # Force: F = -dV/dx = -k * (x - xleft)
                forces[i, 0] = -self.k * (x - self.xleft)
            elif x > self.xright:
                # Right wall: V = k * (x - xright)^2 / 2
                energy += 0.5 * self.k * (x - self.xright)**2
                # Force: F = -dV/dx = -k * (x - xright)
                forces[i, 0] = -self.k * (x - self.xright)
            else:
                # Flat region: V = 0, F = 0
                # energy += 0.0  # Already initialized to 0
                forces[i, 0] = 0.0

        self.results['energy'] = energy
        self.results['forces'] = forces

    def set_parameters(self, **params):
        """Update potential parameters.
        
        Parameters:
            xleft (float): Left boundary
            xright (float): Right boundary  
            k (float): Spring constant
        """
        for key, value in params.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                import warnings
                warnings.warn(f'Ignored unknown parameter "{key}"')
        
        # Validate parameters
        if self.xleft >= self.xright:
            import warnings
            warnings.warn("Setting xleft >= xright in FlatWalls1D potential!")

    def get_boundaries(self):
        """Return the boundaries of the flat region.
        
        Returns:
            tuple: (xleft, xright) boundaries
        """
        return self.xleft, self.xright

    def get_spring_constant(self):
        """Return the spring constant of the walls.
        
        Returns:
            float: Spring constant k
        """
        return self.k

    def potential_at_position(self, x):
        """Calculate potential energy at a single position.
        
        Parameters:
            x (float): Position coordinate
            
        Returns:
            float: Potential energy at position x
        """
        if x < self.xleft:
            return 0.5 * self.k * (x - self.xleft)**2
        elif x > self.xright:
            return 0.5 * self.k * (x - self.xright)**2
        else:
            return 0.0

    def force_at_position(self, x):
        """Calculate force at a single position.
        
        Parameters:
            x (float): Position coordinate
            
        Returns:
            float: Force at position x
        """
        if x < self.xleft:
            return -self.k * (x - self.xleft)
        elif x > self.xright:
            return -self.k * (x - self.xright)
        else:
            return 0.0


class FlatWalls2D(Calculator):
    """ASE-compatible 2D flat walls potential.
    
    Extension of FlatWalls1D to 2D with independent walls in x and y directions.
    
    Parameters:
        xleft, xright (float): x-direction boundaries (default: -0.2, 1.0)
        yleft, yright (float): y-direction boundaries (default: -0.2, 1.0)
        kx, ky (float): Spring constants for x and y walls (default: 100.0, 100.0)
    """
    
    implemented_properties = ['energy', 'forces']

    def __init__(self, xleft=-0.2, xright=1.0, yleft=-0.2, yright=1.0, 
                 kx=100.0, ky=100.0, auto_constrain=True, **kwargs):
        super().__init__(**kwargs)
        self.xleft = xleft
        self.xright = xright
        self.yleft = yleft
        self.yright = yright
        self.kx = kx
        self.ky = ky
        self.auto_constrain = bool(auto_constrain)
        self._constraint_applied = False

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        # Auto-apply constraint (fix z) the first time calculate is called
        if getattr(self, 'auto_constrain', False) and not getattr(self, '_constraint_applied', False):
            atoms.set_constraint(create_2d_constraint(atoms))
            self._constraint_applied = True

        positions = atoms.get_positions()

        energy = 0.0
        forces = np.zeros_like(positions)

        for i in range(len(atoms)):
            x, y = positions[i, 0], positions[i, 1]
            
            # X-direction walls
            if x < self.xleft:
                energy += 0.5 * self.kx * (x - self.xleft)**2
                forces[i, 0] = -self.kx * (x - self.xleft)
            elif x > self.xright:
                energy += 0.5 * self.kx * (x - self.xright)**2
                forces[i, 0] = -self.kx * (x - self.xright)
            else:
                forces[i, 0] = 0.0
            
            # Y-direction walls
            if y < self.yleft:
                energy += 0.5 * self.ky * (y - self.yleft)**2
                forces[i, 1] = -self.ky * (y - self.yleft)
            elif y > self.yright:
                energy += 0.5 * self.ky * (y - self.yright)**2
                forces[i, 1] = -self.ky * (y - self.yright)
            else:
                forces[i, 1] = 0.0

        self.results['energy'] = energy
        self.results['forces'] = forces

    def set_parameters(self, **params):
        """Update potential parameters."""
        for key, value in params.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                import warnings
                warnings.warn(f'Ignored unknown parameter "{key}"')


def create_2d_constraint(atoms):
    """Return a FixCartesian constraint that fixes motion to the x-y plane.

    Usage:
        atoms.set_constraint(create_2d_constraint(atoms))
    """
    # Fix z only (allow x and y motion) for each atom
    return [FixCartesian(i, [False, False, True]) for i in range(len(atoms))]


def create_1d_constraint(atoms):
    """Return a FixCartesian constraint that fixes motion to the x-axis.

    Usage:
        atoms.set_constraint(create_1d_constraint(atoms))
    """
    # Fix y and z (allow x motion only) for each atom
    return [FixCartesian(i, [False, True, True]) for i in range(len(atoms))]
