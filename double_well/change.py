import numpy as np
import os

for i in range(50, 60000):
    path = f"load/{i}/order.txt"
    if not os.path.exists(path):
        break
    data = np.loadtxt(path)
    if data[-1, 1] > 0.5:
        if data[-1, 2] < 0:
            data[:, 2] *= -1
            np.savetxt(path, data)
