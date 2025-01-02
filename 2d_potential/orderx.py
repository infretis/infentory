import numpy as np
from infretis.classes.orderparameter import OrderParameter

class OrderX(OrderParameter):
    def __init__(self, index=0, dim = 0, inter_a=-0.15, inter_b=0.18, energy_a=-9, energy_b=-9):
        super().__init__(description = "OrderX")
        self.inter_a = inter_a
        self.inter_b = inter_b
        self.energy_a = energy_a
        self.energy_b = energy_b
        self.index = 0
        self.dim = 0
        # calculate the potential again because we can't access within ase_engine atm.
        self.g1=1
        self.g2=-10
        self.g3=-10
        self.a1=-30
        self.a2=-3
        self.b1=-30
        self.b2=-3
        self.x0=0.2
        self.y0=0.4

    def vpot(self,x,y):
        V = self.g1 * (x**2 + y**2)**2 \
            + self.g2 * np.exp(self.a1 * (x - self.x0)**2 + self.a2 * (y - self.y0)**2) \
            + self.g3 * np.exp(self.b1 * (x + self.x0)**2 + self.b2 * (y + self.y0)**2)
        return V

    def calculate(self, system):
        pos = system.pos[self.index]
        lmb = pos[self.dim]
        vpot = self.vpot(*pos[:2])

        if lmb < self.inter_a:
            if vpot > self.energy_a:
                lmb = self.inter_a
        elif lmb > self.inter_b:
            if vpot > self.energy_b:
                lmb = self.inter_b

        return [lmb, *pos[:2], vpot]
