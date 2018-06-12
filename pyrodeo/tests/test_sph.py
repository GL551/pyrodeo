#!/usr/bin/python

import numpy as np

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pyrodeo

def test_sph():
    sim = pyrodeo.Simulation.from_geom('sph',
                                       dimensions=[128, 1, 32],
                                       domain=([0.4, 2.5],
                                               [],
                                               [0.5*np.pi-0.15, 0.5*np.pi]),
                                       log_radial=True)

    # Sound speed constant H/r = 0.05
    sim.state.soundspeed = 0.05*sim.state.soundspeed
    sim.param.boundaries[0] = 'reflect'
    sim.param.boundaries[1] = 'periodic'
    sim.param.boundaries[2] = 'reflect'

    # Simulate up to 0.01 orbits
    sim.evolve([0.01*2.0*np.pi], new_file=True)
