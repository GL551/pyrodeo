# -*- coding: utf-8 -*-
"""Simulation parameters used in pyrodeo.

The Param class holds parameters necessary for a pyrodeo simulation.
"""


from __future__ import print_function

import numpy as np
import h5py

class Param(object):
    """Create instance from list of parameters.

    Args:
        param_list ([str, float, float, float, str, str, str]): List of parameters; geometry (string), courant (float), fluxlimiter (float), frame_rotation (float), log_radial (bool), boundaries x (string), boundaries y (string) and boundaries z (string)

    Note:
        The validity of the parameters is not checked.

    The following public attributes are available:

    Attributes:
        geometry (string): 'cart' (Cartesian coordinates), 'sheet' (shearing sheet), 'cyl' (cylindrical coordinates) or 'sph' (spherical coordinates).
        courant (float): Courant number, should be > 0 and < 1.
        limiter_param (float): Limiter parameter. Should be between 1 (minmod) and 2 (superbee).
        min_dens (float): Minimum density to switch to HLL solver to remain positive.
        frame_rotation (float): Frame rotation rate. Ignored in Cartesian coordinates, should be unity in a shearing sheet calculation and corresponds to the angular velocity of the coordinate frame in cylindrical coordinates.
        log_radial (bool): Flag whether to use logarithmic radial coordinates in cylindrical geometry.
        boundaries (string, string, string): Boundary conditions in x y and z: 'nonreflect', 'reflect', or 'periodic'. In shearing sheet mode, the x boundary can be 'shear periodic'.

    """

    def __init__(self, param_list):
        self.geometry       = param_list[0]
        self.courant        = param_list[1]
        self.limiter_param  = param_list[2]
        self.min_dens       = param_list[3]
        self.frame_rotation = param_list[4]
        self.log_radial     = param_list[5]
        self.boundaries     = [param_list[6], param_list[7], param_list[8]]

    @classmethod
    def from_geom(cls,
                  geometry,
                  log_radial=False,
                  boundaries=['reflect', 'reflect', 'reflect']):
        """Initialization from geometry and boundary conditions.

        Construct Parameter object from geometry and boundary conditions. All other parameters are set to standard values. Check if geometry and boundary conditions are valid.

        Args:
            geometry (string): 'cart' (Cartesian coordinates), 'sheet' (shearing sheet), 'cyl' (cylindrical coordinates) or 'sph' (spherical coordinates).
            boundaries (string, string, string): boundary conditions; 'reflect', 'nonreflect' or 'periodic'. In shearing sheet mode, the x boundary can be 'shear periodic'.

        """
        courant = 0.4
        limiter_param = 1.0
        min_dens = 1.0e-6
        frame_rotation = 1.0
        if (geometry != 'cart' and
            geometry != 'sheet' and
            geometry != 'cyl' and
            geometry != 'sph'):
            raise ValueError('Invalid geometry')
        if (log_radial == True and
            geometry != 'cyl' and
            geometry != 'sph'):
            raise ValueError('Can only use logarithmic in cylindrical of spherical coordinates')
        if len(boundaries) != 3:
            raise TypeError('Expected boundaries to have three elements')
        if (boundaries[0] != 'reflect' and
            boundaries[0] != 'nonreflect' and
            boundaries[0] != 'periodic' and
            boundaries[0] != 'shear periodic'):
            raise ValueError('Invalid x boundary')
        if (boundaries[1] != 'reflect' and
            boundaries[1] != 'nonreflect' and
            boundaries[1] != 'periodic'):
            raise ValueError('Invalid y boundary')
        if (boundaries[2] != 'reflect' and
            boundaries[2] != 'nonreflect' and
            boundaries[2] != 'periodic'):
            raise ValueError('Invalid z boundary')

        return cls([geometry,
                    courant,
                    limiter_param,
                    min_dens,
                    frame_rotation,
                    log_radial,
                    boundaries[0],
                    boundaries[1],
                    boundaries[2]])

    def to_list(self):
        """Convert parameters to list.

        Return list of parameters together with list of types. Used for HDF5 output.

        Returns:
            : array of types, list of parameters.

        """
        param_dtype = np.dtype([('geom', h5py.special_dtype(vlen=str)),
                                ('courant', np.float),
                                ('limiter_param', np.float),
                                ('min_dens', np.float),
                                ('frame_rotation', np.float),
                                ('log_radial', np.bool),
                                ('boundary_x', h5py.special_dtype(vlen=str)),
                                ('boundary_y', h5py.special_dtype(vlen=str)),
                                ('boundary_z', h5py.special_dtype(vlen=str))])

        return param_dtype, [(self.geometry,
                              self.courant,
                              self.limiter_param,
                              self.min_dens,
                              self.frame_rotation,
                              self.log_radial,
                              self.boundaries[0],
                              self.boundaries[1],
                              self.boundaries[2])]
