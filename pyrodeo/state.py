# -*- coding: utf-8 -*-
"""Definition of State class used in pyrodeo.

The State class holds density, velocity and sound speed for a pyrodeo simulation
"""

from __future__ import print_function

import numpy as np

class State(object):
    """Construct state holding density, velocity and sound speed for a pyrodeo simulation.

    Args:
        dens (ndarray): 3D ndarray containing density.
        velx (ndarray): 3D ndarray containing x velocity.
        vely (ndarray): 3D ndarray containing y velocity.
        velz (ndarray): 3D ndarray containing z velocity.
        soundspeed (ndarray): 3D ndarray containing sound speed.

    Note:
        No checks are performed whether density, velocity and sound speed are valid arrays. They should all have the same shape, the same as the arrays of :class:`.Coordinates`.

    The following public attributes are available:

    Attributes:
        dens (ndarray): 3D ndarray containing density.
        velx (ndarray): 3D ndarray containing x velocity.
        vely (ndarray): 3D ndarray containing y velocity.
        velz (ndarray): 3D ndarray containing z velocity.
        soundspeed (ndarray): 3D ndarray containing sound speed.
        no_ghost (ndarray): 3D ndarray flagging whether a cell is a ghost cell (=0) or an internal cell (=1)

    """

    def __init__(self, dens, velx, vely, velz, soundspeed):
        self.dens = dens
        self.velx = velx
        self.vely = vely
        self.velz = velz
        self.soundspeed = soundspeed

        self.no_ghost = np.full(np.shape(dens), 1)
        if len(dens[:,0,0]) > 1:
            self.no_ghost[:2,:,:]  = 0
            self.no_ghost[-2:,:,:] = 0
        if len(dens[0,:,0]) > 1:
            self.no_ghost[:,:2,:]  = 0
            self.no_ghost[:,-2:,:] = 0
        if len(dens[0,0,:]) > 1:
            self.no_ghost[:,:,:2]  = 0
            self.no_ghost[:,:,-2:] = 0

    @classmethod
    def from_dims(cls, dims):
        """Construct State from grid dimensions.

        Construct State given grid dimensions, creating arrays of the correct size with standard (physical) values.

        Args:
            dims (int, int, int): Dimensions of the grid in x, y and z.

        """
        if len(dims) != 3:
            raise TypeError('Expexted dimensions to have two elements')
        if (dims[0] < 1 or dims[1] < 1 or dims[2] < 1):
            raise ValueError('Need all dimensions to be larger than zero')

        dens = np.full(dims, 1.0)
        velx = np.full(dims, 0.0)
        vely = np.full(dims, 0.0)
        velz = np.full(dims, 0.0)
        soundspeed = np.full(dims, 1.0)

        return cls(dens, velx, vely, velz, soundspeed)

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
        velz = np.copy(other_state.velz)
        soundspeed = np.copy(other_state.soundspeed)

        return cls(dens, velx, vely, velz, soundspeed)

    def transpose(self, axis_order):
        """Change axis order for all fields."""
        self.dens = np.transpose(self.dens, axis_order)
        self.velx = np.transpose(self.velx, axis_order)
        self.vely = np.transpose(self.vely, axis_order)
        self.velz = np.transpose(self.velz, axis_order)
        self.soundspeed = np.transpose(self.soundspeed, axis_order)
        self.no_ghost = np.transpose(self.no_ghost, axis_order)

    def swap_velocities(self, dim):
        """Swap two velocities."""
        if dim == 1:
            tmp = self.velx
            self.velx = self.vely
            self.vely = tmp
        if dim == 2:
            tmp = self.velx
            self.velx = self.velz
            self.velz = tmp
