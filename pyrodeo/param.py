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
        param_list ([str, float, float, float, str, str, str]): List of parameters; geometry (string), courant (float), fluxlimiter (float), frame_rotation (float), log_radial (bool), boundaries x in (string), boundaries x out (string), boundaries y in (string), boundaries y out (string), boundaries z in (string) and boundaries z out (string)

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
        boundaries ((str,str), (str,str), (str,str)): Boundary conditions (in and out) in x y and z: 'nonreflecting', 'closed', 'symmetric', or 'periodic'. In shearing sheet mode, the x boundary can be 'shear periodic'.

    """

    def __init__(self, param_list):
        self.geometry       = param_list[0]
        self.courant        = param_list[1]
        self.limiter_param  = param_list[2]
        self.min_dens       = param_list[3]
        self.frame_rotation = param_list[4]
        self.log_radial     = param_list[5]
        self.boundaries     = [[param_list[6], param_list[7]],
                               [param_list[8], param_list[9]],
                               [param_list[10], param_list[11]]]

    @classmethod
    def from_geom(cls,
                  geometry,
                  log_radial=False,
                  boundaries=[['closed','closed'],
                              ['closed','closed'],
                              ['closed','closed']]):
        """Initialization from geometry and boundary conditions.

        Construct Parameter object from geometry and boundary conditions. All other parameters are set to standard values. Check if geometry and boundary conditions are valid.

        Args:
            geometry (string): 'cart' (Cartesian coordinates), 'sheet' (shearing sheet), 'cyl' (cylindrical coordinates) or 'sph' (spherical coordinates).
            boundaries ((str,str), (str,str), (str,str)): boundary conditions; 'closed', 'nonreflecting', 'symmetric', or 'periodic'. In shearing sheet mode, the x boundary can be 'shear periodic'.

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
        for bc in boundaries[0]:
            if (bc != 'closed' and
                bc != 'nonreflecting' and
                bc != 'periodic' and
                bc != 'shear periodic' and
                bc != 'symmetric'):
                raise ValueError('Invalid x boundary')
        for bc in boundaries[1]:
            if (bc != 'closed' and
                bc != 'nonreflecting' and
                bc != 'periodic' and
                bc != 'symmetric'):
                raise ValueError('Invalid y boundary')
        for bc in boundaries[2]:
            if (bc != 'closed' and
                bc != 'nonreflecting' and
                bc != 'periodic' and
                bc != 'symmetric'):
                raise ValueError('Invalid z boundary')
        if (boundaries[0][0] == 'shear periodic' and
            geometry != 'sheet'):
            raise ValueError('Can only use shear periodic boundaries in shearing sheet geometry')

        return cls([geometry,
                    courant,
                    limiter_param,
                    min_dens,
                    frame_rotation,
                    log_radial,
                    boundaries[0][0],
                    boundaries[0][1],
                    boundaries[1][0],
                    boundaries[1][1],
                    boundaries[2][0],
                    boundaries[2][1]])

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
                                ('boundary_x_in', h5py.special_dtype(vlen=str)),
                                ('boundary_x_out', h5py.special_dtype(vlen=str)),
                                ('boundary_y_in', h5py.special_dtype(vlen=str)),
                                ('boundary_y_out', h5py.special_dtype(vlen=str)),
                                ('boundary_z_in', h5py.special_dtype(vlen=str)),
                                ('boundary_z_out', h5py.special_dtype(vlen=str))])

        return param_dtype, [(self.geometry,
                              self.courant,
                              self.limiter_param,
                              self.min_dens,
                              self.frame_rotation,
                              self.log_radial,
                              self.boundaries[0][0],
                              self.boundaries[0][1],
                              self.boundaries[1][0],
                              self.boundaries[1][1],
                              self.boundaries[2][0],
                              self.boundaries[2][1])]
