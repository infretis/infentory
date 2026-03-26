import numpy as np
from infretis.classes.orderparameter import OrderParameter

class OrderX(OrderParameter):
    def __init__(self, a, b, c):
        super().__init__(description = "OrderX")
        self.a = a
        self.b = b
        self.c = c
        self.angles = np.deg2rad(np.arange(0, 91, 5))

    def calculate(self, system):
        pos = system.pos[0]
        x, y = pos[0], pos[1]
        proj = x * np.cos(self.angles) + y * np.sin(self.angles)

        return [x, y] + list(proj)
