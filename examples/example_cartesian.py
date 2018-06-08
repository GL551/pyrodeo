#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pyrodeo

# Create simulation with default resolution and domain
sim = pyrodeo.Simulation.from_geom('cart')

# The basic state will have density and sound speed unity everywhere,
# and the velocity will be zero. In order to create a simple shock tube,
# now set density to 1/10 for x > 0.
sel = np.where(sim.coords.x > 0.0)
sim.state.dens[sel] = 0.1

# Evolve until t = 0.25
sim.evolve([0.25], new_file=True)

# Plot results
plt.plot(sim.coords.x, sim.state.dens)
plt.show()
