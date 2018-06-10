#!/usr/bin/python

import numpy as np
import pyrodeo

# Extra source terms: planet gravity
def planet_source(t, dt, coords, state, planetParam):
    # Mass ratio planet/star
    mp = planetParam[0]
    # Softening length planet potential
    eps = planetParam[1]

    # Coordinates
    r = coords.x
    p = coords.y

    # Planet coordinates
    rp = 1.0
    pp = 0.0

    # Distance to the planet
    dist = np.sqrt(r*r + rp*rp - 2.0*r*rp*np.cos(p - pp) + eps*eps)

    # Potential gradient
    dpotdr = mp*(r - rp*np.cos(p - pp))/(dist*dist*dist)
    dpotdp = mp*r*rp*np.sin(p - pp)/(dist*dist*dist)

    # Indirect term
    dpotdr += mp*np.cos(p - pp)/(rp*rp)
    dpotdp -= mp*r*np.sin(p - pp)/(rp*rp)

    # Resulting source term
    source_velx = -dpotdr
    source_vely = -dpotdp/(r*r)

    # Damping boundary conditions
    Rin = 100.0*(r - 0.5)*(r - 0.5)
    Rin[np.where(r > 0.5)] = 0.0
    Rout = (r - 2.1)*(r - 2.1)/(0.4*0.4)
    Rout[np.where(r < 2.1)] = 0.0
    R = (Rin + Rout)*np.power(r, -1.5)

    # Damp towards initial state
    source_dens = -(state.dens - 1.0)*R
    source_velx -= state.velx*R
    source_vely -= state.vely*R

    # Integrate extra source terms
    state.dens += dt*source_dens*state.no_ghost
    state.velx += dt*source_velx*state.no_ghost
    state.vely += dt*source_vely*state.no_ghost

sim = pyrodeo.Simulation.from_geom('cyl',
                                   dimensions=[128, 384],
                                   domain=([0.4, 2.5], [-np.pi, np.pi]))

# Sound speed constant H/r = 0.05
sim.state.soundspeed = 0.05*sim.state.soundspeed/np.sqrt(sim.coords.x)
sim.param.boundaries[0] = 'reflect'
sim.param.boundaries[1] = 'periodic'

# Simulate a Jupiter planet up to 100 orbits
sim.evolve(2.0*np.pi*np.array([1.0,2.0,5.0,10.0,20.0,50.0,100.0]),
           planet_source, (0.001, 0.6*0.05), new_file=True)
