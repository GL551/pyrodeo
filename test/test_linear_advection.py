#!/usr/bin/python

import numpy as np

from pyrodeo.state import State
from pyrodeo.coords import Coordinates
import pyrodeo.linear_advection as pr

def linear_advection_error():
    coords = Coordinates.from_dims(dimensions=(100, 1),
                                   domain=([-0.5, 0.5], []))
    state = State.from_dims((104,1))

    state.dens += np.cos(2.0*np.pi*coords.x)
    advection_velocity = 0.0*coords.x + 1.0

    la = pr.LinearAdvection(advection_velocity, 1.0)

    maxn = 17
    for n in range(maxn):
        la.step(1.0/maxn, coords.dxy[0], state)

    dens_error = np.sum(np.abs(state.dens - 1.0 -
                               np.cos(2.0*np.pi*coords.x)))

    return dens_error

def test_linear_advection():
    assert linear_advection_error() < 0.02
