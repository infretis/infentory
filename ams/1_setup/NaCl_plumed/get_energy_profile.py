#!/usr/bin/env amspython
from scm.plams import *
import os
import numpy as np
import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import metadynminer
import sys

## input

N_interfaces = int(sys.argv[1])
dE_cap = 1
e_shift_lambda_a = 4.5
cv_shift_lambda_b = .05
cv_shift_cap = 0.0

def derivative(x,y):
    shape = y.shape[0]
    x = x[0:shape]
    dy=np.diff(y,1)
    dx=np.diff(x,1)
    y_der=dy/dx
    sign = np.sign(y_der)
    signchange = ((np.roll(sign, 1) - sign) != 0)
    signchange = signchange[1:]
    # signchange = np.append(signchange, False)
    return y_der, signchange

def interface_positions(e_interface, slope_cv, slope_e):
    positions = []
    e_interface = list(e_interface)

    print(slope_cv)
    for i in e_interface:
        position = slope_cv[slope_e.searchsorted(i, 'left')]
        positions.append(position)
    return positions

hillsfile = metadynminer.Hills(name="HILLS", periodic=[False])
print(hillsfile.hills)
fes = metadynminer.Fes(hills=hillsfile, resolution=1024, original=True)#, cv1range=[0, 90])
print(fes)
# fes.plot(png_name="fes.png")
print(fes.fes.shape)
print(fes.cv1)
profile  = fes.fes
cv_range = fes.cv1max-fes.cv1min
cv1min = fes.cv1min - cv_range*0.15
cv1max = fes.cv1max + cv_range*0.15
print(fes.cv1)
print(fes.cv1max, fes.cv1min)
cv = np.linspace(cv1min, cv1max, fes.res)


fig, ax = plt.subplots()
ax.plot(cv, fes.fes)
try:
    cv_der = cv[1:-1]
    profile_short = profile[1:-1]
    prof_1der, pos_extrema = derivative(cv, profile)
    prof_2der, _ = derivative(cv, prof_1der)
    type_extrema = prof_2der[pos_extrema]
    fe_extrema = profile_short[pos_extrema]
    pos_extrema = cv_der[pos_extrema]
    print(pos_extrema, type_extrema, fe_extrema)
    barrier = type_extrema < 0
    state_A = type_extrema[0] 
    pos_state_a = pos_extrema[0]
    pos_barrier = pos_extrema[1]
    pos_state_b = pos_extrema[2]
    lambda_cap = pos_barrier - cv_range * cv_shift_cap #+ cv_range * 0.05
    e_lambda_i = fe_extrema[1] - dE_cap
    e_lambda_a = fe_extrema[0] + e_shift_lambda_a
    lambda_b = pos_state_b - cv_range * cv_shift_lambda_b

    e_steps = np.linspace(e_lambda_a, e_lambda_i, N_interfaces)
    dE = np.round(e_steps[1] - e_steps[0], 2)
    slope = np.where(np.logical_and(cv >= pos_state_a, cv <= pos_barrier))
    slope_cv = cv[slope]
    slope_e = profile[slope]
    interfaces = interface_positions(e_steps, slope_cv, slope_e)
    # lambda_b = - interfaces[0]

    interfaces.append(lambda_b)
    print(interfaces)
    print(e_lambda_a)
    print(pos_state_a, pos_state_b, pos_barrier, lambda_cap)

    for i, x in enumerate(interfaces):
        if i == 0:
            plt.axvline(x=x, color='blue', label='lambda_a')
        elif i == len(interfaces)-1:
            plt.axvline(x=x, color='red', label='lambda_b')
        else:
            plt.axvline(x=x, color='grey', linestyle='dotted')

    plt.axvline(x=lambda_cap, color='black', label='cap')
    plt.legend()
    plt.xlabel("CV")
    plt.ylabel("Free Energy [kJ/mol]")
    plt.text(pos_state_a, fe_extrema[0], "A", verticalalignment='bottom', 
            horizontalalignment='center', fontsize=14 )
    plt.text(pos_state_b, fe_extrema[2], "B", verticalalignment='bottom',
            horizontalalignment='center', fontsize=14 )
    plt.text(0.99, 0.345, r'$\Delta F_{int}=$'+f'{dE} kj/mol', 
            horizontalalignment='right', verticalalignment='top', 
            transform = ax.transAxes, fontsize=11)
    plt.text(0.99, 0.3, r'$\Delta F_{cap}=$'+f' {dE_cap} kj/mol', 
            horizontalalignment='right', verticalalignment='top', 
            transform = ax.transAxes, fontsize=11)
    plt.text(0.99, 0.255, r'$N_{int}=$'+f'{N_interfaces}', 
            horizontalalignment='right', verticalalignment='top', 
            transform = ax.transAxes, fontsize=11)

    moves = ['sh', 'sh']
    for i in range(len(interfaces)-2):
        moves.append('wf')

    interfaces = [ round(elem, 2)  for elem in interfaces]
    lambda_cap = round(lambda_cap, 2)

    with open(f'interfaces_{N_interfaces}.toml', 'w') as f:
        f.writelines('[simulation]\n')
        f.writelines(f'interfaces = {interfaces}\n')
        f.writelines(f'shooting_moves = {moves}\n\n')
        f.writelines('[simulation.tis_set]\n')
        f.writelines(f'interface_cap = {lambda_cap}')

    minima = metadynminer.Minima(fes)
    plt.savefig(f'profile_{N_interfaces}.png')
    plt.show()
except Exception:
    plt.show()
