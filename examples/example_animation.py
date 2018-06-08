#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import h5py

# Return density at checkpoint n
def density_at_checkpoint(file_name, n):
    dens = None
    with h5py.File(file_name, 'r') as hf:
        maxn = len(hf.keys()) - 2
        if n >= maxn:
            n = maxn - 1

        s = "checkpoint{}".format(n)
        g = hf.get(s)
        dens = np.array(g.get('dens'))

    return dens

fig = plt.figure()

# Show initial conditions (checkpoint 0)
file_name = './rodeo.h5'
dens = density_at_checkpoint(file_name, 0)
im = plt.imshow(dens, animated=True)

# Get next checkpoint
n = 0
def updatefig(*args):
    global n

    dens = density_at_checkpoint(file_name, n)
    im.set_array(dens)
    n += 1

    return im,

# Animate!
ani = animation.FuncAnimation(fig, updatefig, interval=50, blit=True)
plt.show()
