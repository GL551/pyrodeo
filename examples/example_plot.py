#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import h5py

with h5py.File('./rodeo.h5', "r") as hf:
    # Select last available checkpoint
    last_checkpoint = None
    for k in hf.keys():
        if (k != 'coords' and k != 'param'):
            last_checkpoint = k

    # Get x coordinate
    gc = hf.get('coords')
    x = np.array(gc.get('x'))

    # Get density
    g = hf.get(last_checkpoint)
    dens = np.array(g.get('dens'))
    # Simulation time at checkpoint
    t = g.attrs['time']

    print('Plotting ' + last_checkpoint + ' at t = {}'.format(t))

    plt.plot(x[:,0], np.mean(dens, axis=1))

    plt.show()
