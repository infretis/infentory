import numpy as np
from ase import units
from ase.calculators.calculator import Calculator, all_changes

class Potential2D(Calculator):
    """
    The 2D potential from:
        https://pyretis.org/current/examples/examples-2d-hysteresis.html

    For unit conversion:
        https://www.pyretis.org/current/user/units.html
        https://docs.lammps.org/units.html

    """
    def __init__(self, g1=1, g2=-10, g3=-10, a1=-30, a2=-3, b1=-30, b2=-3, x0=0.2, y0=0.4):
        self.implemented_properties = ["energy", "forces"]
        super().__init__()
        self.g1 = g1
        self.g2 = g2
        self.g3 = g3
        self.a1 = a1
        self.a2 = a2
        self.b1 = b1
        self.b2 = b2
        self.x0 = x0
        self.y0 = y0

    def calculate(self, atoms = None, system_changes = all_changes, properties = ["energy", "forces"]):
        x, y = atoms.positions[0,:2]

        V = self.g1 * (x**2 + y**2)**2 \
            + self.g2 * np.exp(self.a1 * (x - self.x0)**2 + self.a2 * (y - self.y0)**2) \
            + self.g3 * np.exp(self.b1 * (x + self.x0)**2 + self.b2 * (y + self.y0)**2)

        dV_dx = 4 * self.g1 * x * (x**2 + y**2) \
                + 2 * self.g2 * self.a1 * (x - self.x0) * np.exp(self.a1 * (x - self.x0)**2 + self.a2 * (y - self.y0)**2) \
                + 2 * self.g3 * self.b1 * (x + self.x0) * np.exp(self.b1 * (x + self.x0)**2 + self.b2 * (y + self.y0)**2)

        dV_dy = 4 * self.g1 * y * (x**2 + y**2) \
                + 2 * self.g2 * self.a2 * (y - self.y0) * np.exp(self.a1 * (x - self.x0)**2 + self.a2 * (y - self.y0)**2) \
                + 2 * self.g3 * self.b2 * (y + self.y0) * np.exp(self.b1 * (x + self.x0)**2 + self.b2 * (y + self.y0)**2)

        F_x = -dV_dx
        F_y = -dV_dy

        self.results = {}
        self.results["energy"] = V
        self.results["forces"] = np.array([[F_x, F_y, 0.0]])

def plotpotential():
    import matplotlib.pyplot as plt
    from ase.atoms import Atoms
    from potential import Potential2D

    a = Atoms("H")
    a.calc = Potential2D()
    N = 75
    x0 = np.linspace(-0.5,0.5,N);y0 = np.linspace(-1,1,N)
    X = np.meshgrid(x0,y0)
    V = []
    for i,x in enumerate(X[0]):
        for j,y in enumerate(X[1]):
            a.positions[0,0] = x[i]
            a.positions[0,1] = y[i]
            V.append(a.get_potential_energy())

    plt.scatter(X[0].flatten(),X[1].flatten(), c = V)
    plt.show()
