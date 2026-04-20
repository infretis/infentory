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
    def __init__(self, a, b, c):
        self.implemented_properties = ["energy", "forces"]
        super().__init__()
        self.a = a
        self.b = b
        self.c = c

    def calculate(self, atoms = None, system_changes = all_changes, properties = ["energy", "forces"]):
        x, y = atoms.positions[0,:2]

        v_pot = self.a*(y**4 + x**4 -self.b*x**2 - self.c*x)

        dv_dx = self.a*(4*x**3 - 2*self.b*x - self.c)  # x derivative
        dv_dy = 4*self.a*y**3                    # y derivative

        self.results = {}
        self.results["energy"] = v_pot
        self.results["forces"] = np.array([[-dv_dx, -dv_dy, 0.0]])

def plotpotential():
    import matplotlib.pyplot as plt
    from ase.atoms import Atoms
    from potential import Potential2D

    a = Atoms("H")
    a.calc = Potential2D(1.3128, 2.0, 0.25)
    N = 75
    x0 = np.linspace(-1.5,1.5,N);y0 = np.linspace(-1.5,1.5,N)
    X = np.meshgrid(x0,y0)
    V = np.zeros((N, N))
    for i,x in enumerate(X[0]):
        for j,y in enumerate(X[1]):
            a.positions[0,0] = x[i]
            a.positions[0,1] = y[i]
            V[i, j] = a.get_potential_energy()

    plt.contour(x0, y0, V.T, levels=20)
    plt.show()

# plotpotential()
