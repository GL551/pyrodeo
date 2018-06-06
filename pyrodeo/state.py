# -*- coding: utf-8 -*-
"""Definition of State class used in pyrodeo.

The State class holds density, velocity and sound speed for a pyrodeo simulation
"""

from __future__ import print_function

import numpy as np

class State(object):
    """Construct state holding density, velocity and sound speed for a pyrodeo simulation.

    Args:
        dens (ndarray): 2D ndarray containing density.
        velx (ndarray): 2D ndarray containing x velocity.
        vely (ndarray): 2D ndarray containing y velocity.
        soundspeed (ndarray): 2D ndarray containing sound speed.

    Note:
        No checks are performed whether density, velocity and sound speed are valid arrays. They should all have the same shape, the same as the arrays of :class:`.Coordinates`.

    The following public attributes are available:

    Attributes:
        dens (ndarray): 2D ndarray containing density.
        velx (ndarray): 2D ndarray containing x velocity.
        vely (ndarray): 2D ndarray containing y velocity.
        soundspeed (ndarray): 2D ndarray containing sound speed.

    """

    def __init__(self, dens, velx, vely, soundspeed):
        self.dens = dens
        self.velx = velx
        self.vely = vely
        self.soundspeed = soundspeed

    @classmethod
    def from_dims(cls, dims):
        """Construct State from grid dimensions.

        Construct State given grid dimensions, creating arrays of the correct size with standard (physical) values.

        Args:
            dims (int, int): Dimensions of the grid in x and y

        """
        if len(dims) != 2:
            raise TypeError('Expexted dimensions to have two elements')
        if (dims[0] < 1 or dims[1] < 1):
            raise ValueError('Need both dimensions to be larger than zero')

        dens = np.full(dims, 1.0)
        velx = np.full(dims, 0.0)
        vely = np.full(dims, 0.0)
        soundspeed = np.full(dims, 1.0)

        return cls(dens, velx, vely, soundspeed)

    @classmethod
    def copy(cls, other_state):
        """Construct state from other State.

        Set this instance of State equal to an other State, performing an explicit copy.

        Args:
            other_state (State): State from which to copy.

        """
        dens = np.copy(other_state.dens)
        velx = np.copy(other_state.velx)
        vely = np.copy(other_state.vely)
        soundspeed = np.copy(other_state.soundspeed)

        return cls(dens, velx, vely, soundspeed)

    def transpose(self):
        """Transpose all fields."""
        self.dens = np.transpose(self.dens)
        self.velx = np.transpose(self.velx)
        self.vely = np.transpose(self.vely)
        self.soundspeed = np.transpose(self.soundspeed)

    def swap_velocities(self):
        """Swap x and y velocities."""
        tmp = self.velx
        self.velx = self.vely
        self.vely = tmp
