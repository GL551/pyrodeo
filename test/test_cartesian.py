#!/usr/bin/python

import numpy as np

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pyrodeo.simulation as prs

def test_cartesian():
    # Create simulation with default resolution and domain
    sim = prs.Simulation.from_geom('cart')

    # The basic state will have density and sound speed unity everywhere,
    # and the velocity will be zero. In order to create a simple shock tube,
    # now set density to 1/10 for x > 0.
    sel = np.where(sim.coords.x > 0.0)
    sim.state.dens[sel] = 0.1

    # Evolve until t = 0.25
    sim.evolve([0.25], new_file=True)
