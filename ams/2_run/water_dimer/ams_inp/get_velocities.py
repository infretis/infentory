from scm.plams import RKFTrajectoryFile, KFFile, Units
import numpy as np

rkf = RKFTrajectoryFile('ams.results/ams.rkf', mode='rb')
rkf.store_mddata()
rkf.store_historydata()
mol = rkf.get_plamsmol()
# kf = KFFile('ams.rkf')

# print(kf.read())
velocities = []
print(rkf.get_length())
for i in range(rkf.get_length()):
    if i > 0:
        print(i)
        rkf.read_frame(i, molecule=mol)   
        velocities.append(Units.convert(rkf.mddata['Velocities'], 'bohr', 'nm'))

print(np.array(velocities).shape)
vels = np.array(velocities).flatten() * 1000 
# make shape (N, 3)
# vels = vels.reshape(-1, 3)
# print(vels.shape)
# # get rowwise pythagorean
# vels = np.linalg.norm(vels, axis=1)
# self.dist_unit = "nm"
#         self.time_unit = "ps"
# get standard deviation
print(vels)
print(vels.shape, rkf.get_length())
std = np.std(vels)
print(std, np.average(vels))

# plot histogram of velocities
import matplotlib.pyplot as plt
plt.hist(vels, bins=100)
plt.savefig('velocities.png')
plt.show()


