# -*- coding: utf-8 -*-
"""Defines the Simulation class holding a pyrodeo simulation
"""

from __future__ import print_function

import numpy as np
import h5py
from pyrodeo.state import State
from pyrodeo.coords import Coordinates
from pyrodeo.param import Param
from pyrodeo.hydro import Hydro

class Simulation(object):
    """Construct pyrodeo simulation.

    Construct simulation from existing instances of :class:`.Param`, :class:`.Coordinates` and :class:`.State` and given simulation time.

    Args:
        param (:class:`.Param`): Valid Param object, containing simulation parameters.
        coords (:class:`.Coordinates`): Valid Coordinates object, containing x and y coordinates
        state (:class:`.State`): Valid State object, containing current density, velocity and sound speed.
        t (float): Simulation time.
        direc (:obj:`str`, optional): Output directory.

    Note:
        While it is possible to set up a simulation using the basic constructor, no checks are performed whether the instances of :class:`.Param`, :class:`.Coordinates` and :class:`.State` are valid. It is safer to use :meth:`.from_geom` which will set up the necessary valid instances.

    The following attributes and methods are available:

    Attributes:
        state (:class:`.State`): State object holding density, velocity and sound speed.
        coords (:class:`.Coordinates`): Coordinates object holding x and y coordinates.
        param (:class:`.Param`): Param object holding simulation parameters.
        t (float): Current simulation time.
        direc (string): Output directory.

    """

    def __init__(self, param, coords, state, t, direc='./'):
        self.state = State(state.dens, state.velx, state.vely, state.soundspeed)
        self.coords = Coordinates(coords.x, coords.y)
        self.param = Param(param)
        self.t = t
        self.direc = direc

    @classmethod
    def from_geom(cls,
                  geometry,
                  dimensions=(100, 1),
                  domain=([-0.5, 0.5], []),
                  direc='./'):
        """Construct simulation from geometry.

        Construct simulation class from geometry, grid dimensions and domain size. All other parameters will be set to defaults. Simulation time is set to zero.

        Args:
            geometry (string): 'cart' (Cartesian coordinates), 'sheet' (shearing sheet), or 'cyl' (cylindrical coordinates).
            dimensions (:obj:`(int,int)`,optional): Grid dimensions in x and y.
            domain (:obj:`[(float,float),(float,float)]`,optional): Domain boundaries in x and y.
            direc (:obj:`str`,optional): Output directory, defaults to current directory.

        """
        dt, param = Param.from_geom(geometry).to_list()

        coords = Coordinates.from_dims(dimensions, domain)

        state = State.from_dims((len(coords.x[:,0]),
                                 len(coords.x[0,:])))
        t = 0.0

        return cls(param[0], coords, state, t, direc=direc)

    @classmethod
    def from_checkpoint(cls, direc, n=None):
        """Constructor for simulation class from checkpoint.

        Construct simulation class from previously saved checkpoint.

        Args:
            direc (str): Path to 'rodeo.h5'. This will be the output directory as well.
            n (:obj:`int`,optional): Index of checkpoint to restore. If left None, restore last saved checkpoint.
        """
        with h5py.File(direc + '/rodeo.h5', "r") as hf:
            # Restart from last checkpoint
            if n is None:
                g = hf.keys()
                n = len(g) - 3

            gc = hf.get("coords")
            x = np.array(gc.get('x'))
            y = np.array(gc.get('y'))

            param = hf.get("param")

            g = hf.get("checkpoint{}".format(n))
            dens = np.array(g.get('dens'))
            velx = np.array(g.get('velx'))
            vely = np.array(g.get('vely'))
            soundspeed = np.array(g.get('soundspeed'))

            coords = Coordinates(x, y)
            state = State(dens, velx, vely, soundspeed)

            t = g.attrs['time']
            print("Restoring from checkpoint at t = {}".format(t))

            return cls(param[0], coords, state, t, direc)

    def __str__(self):
        """Information on simulation."""
        return 'Simulation class, part of pyrodeo.\n \
    Current time: {}\n'.format(self.t) + str(self.coords) + str(self.param)

    def evolve(self, checkpoints,
               source_func = None, source_param = None,
               new_file=False):
        """Evolve simulation over a list of checkpoints.

        Args:
            checkpoints (ndarray): List of times to save checkpoints.
            Source (:obj:`callable`, optional): Extra source terms. Must be of the form f(t, coords, state, source_param) and return three source terms of the same shape as state.dens: density source term, x velocity source term and y velocity source term.
            source_param (:obj:`ndarray`, optional): Parameters for extra source term. Will be passed to the Source function.
            new_file (:obj:`bool`,optional): If true, create new output file 'rodeo.h5', otherwise append to file if it exists.

        """
        solver = Hydro(self.param, self.coords)
        if new_file is True:
            self.checkpoint(new_file=True)
        for tmax in checkpoints:
            if self.t < tmax:
                self.t, self.state = \
                    solver.evolve(self.t,
                                  tmax,
                                  self.coords,
                                  self.param,
                                  self.state,
                                  source_func,
                                  source_param)
                self.checkpoint()
            else:
                print("Skipping checkpoint t = {}".format(tmax))

    def checkpoint(self, new_file=False):
        """Save checkpoint in 'rodeo.h5'

        Args:
            new_file (:obj:`bool`,optional): If true, create new output file 'rodeo.h5', otherwise append to file if it exists.

        """
        mode = "a"
        if new_file is True:
            mode = "w"

        print("Checkpoint at t = {}".format(self.t))

        with h5py.File(self.direc + '/rodeo.h5', mode) as hf:
            n = 0

            # Write coordinates and parameters
            if new_file is True:
                gc = hf.create_group("coords")
                gc.create_dataset("x", data=self.coords.x)
                gc.create_dataset("y", data=self.coords.y)

                param_dtype, param = self.param.to_list()
                param_array = np.array(param, dtype=param_dtype)
                param_dset = \
                  hf.create_dataset("param", (len(param),), dtype=param_dtype)
                param_dset[...] = param_array
            else:
                # Number of previous saves + 1
                g = hf.keys()
                n = len(g) - 2

            # Write current data
            g1 = hf.create_group("checkpoint{}".format(n))
            g1.attrs['time'] = self.t
            g1.create_dataset("dens", data=self.state.dens)
            g1.create_dataset('velx', data=self.state.velx)
            g1.create_dataset('vely', data=self.state.vely)
            g1.create_dataset('soundspeed', data=self.state.soundspeed)
