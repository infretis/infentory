import glob
import os
import shutil

import numpy as np
import matplotlib.pyplot as plt

import MDAnalysis as mda
import tomli


def gen_init_paths():
    with open("../infretis.toml", "rb") as toml_file:
        tdict= tomli.load(toml_file)
    intfs = tdict["simulation"]["interfaces"]

    orderp = np.loadtxt('order.txt')
    traj = mda.Universe('traj_comp.xtc')
    os.makedirs('load_init', exist_ok=True)
    lengths = [5*(i+1) for i in range(len(intfs))]

    print(intfs)
    for idx, intf in enumerate(intfs[:-1]):
        # aib = [intfs[0], intf, intfs[-1]]
        # print('0', aib)
        aib = [intfs[0] if idx == 0 else intfs[idx-1], intf, intfs[idx+1]]
        # print('1', aib0)
        folder = f'load_init/{idx+1}'
        os.makedirs(folder, exist_ok=True)
        os.makedirs(folder+'/accepted', exist_ok=True)
        gen_paths(aib, orderp, traj, folder, plen=lengths[idx])

    # gen for minus
    aib = [intfs[0], intfs[0], intfs[0]]
    folder = 'load_init/0'
    os.makedirs(folder, exist_ok=True)
    os.makedirs(folder+'/accepted', exist_ok=True)
    gen_paths(aib, orderp, traj, folder, minus=True)


def gen_paths(intfs_aib, orderp, traj, folder, minus=False, plen=10):
    # print(orderp[:, 1])
    # print(orderp[:, 1])
    print('cheese')
    argsort = np.argsort(orderp[:, 1])
    und_a = orderp[:, 1] < intfs_aib[0]
    ove_a = orderp[:, 1] > intfs_aib[0]
    bet_ib = np.logical_and(intfs_aib[1] < orderp[:, 1],
                            orderp[:, 1] < intfs_aib[2])

    print('dow', intfs_aib)
    # idxes
    below = np.where(und_a == 1)[0]
    inbet = np.where(bet_ib == 1)[0]
    above = np.where(ove_a == 1)[0]

    below_args = np.argsort(orderp[:, 1][below])
    inbet_args = np.argsort(orderp[:, 1][inbet])
    above_args = np.argsort(orderp[:, 1][above])
    below_s = below[below_args]
    inbet_s = inbet[inbet_args][:plen][::1]
    print('zum', len(inbet_s))
    above_s = above[above_args]

    if minus:
        traj_idxes = [above_s[0]] + list(below_s[-10:]) + [above_s[0]]
    else:
        print('crown', below_s)
        traj_idxes = [below_s[-1]] + list(inbet_s) + [below_s[-1]]

    trajfile = os.path.join(folder, 'accepted', "traj.trr")
    orderfile = os.path.join(folder, "order.txt")
    trajtxtfile = os.path.join(folder, "traj.txt")

    N = len(traj_idxes)
    np.savetxt(
        orderfile,
        np.c_[orderp[:N, 0], orderp[traj_idxes, 1]],
        header=f"{'time':>10} {'order':>15}",
        fmt=["%10.d", "%15.4f"],
    )
    np.savetxt(
        trajtxtfile,
        np.c_[
            [str(i) for i in range(N)],
            ["traj.trr" for i in range(N)],
            [str(i) for i in range(N)],
            [str(1) for i in range(N)],
        ],
        header=f"{'time':>10} {'trajfile':>15} {'index':>10} {'vel':>5}",
        fmt=["%10s", "%15s", "%10s", "%5s"],
    )
    all_at = traj.select_atoms("all")
    with mda.Writer(trajfile, all_at.n_atoms) as W:
        for ts in traj.trajectory[traj_idxes]:
            W.write(all_at)
    # traj 
    # atoms = [i.name for i in traj.atoms]
    # with open(trajfile, 'w') as write:
    #     for idx in traj_idxes:
    #         write.write(f'{len(atoms)}\n')
    #         write.write('whada\n')
    #         pos = traj.trajectory[idx]
    #         for atom, xyz in zip(atoms, pos):
    #             xyz_s = [atom] + [f'{i:.8f}' for i in xyz]
    #             write.write('\t'.join(xyz_s) + '\n')


    print('cow', len(traj_idxes))

    print(sum(bet_ib))
    print(bet_ib)


gen_init_paths()
