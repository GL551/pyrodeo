# -*- coding: utf-8 -*-
"""Coordinate class used in pyrodeo.

The Coordinate class hold the x, y and z coordinates as 3D ndarrays. In addition, it holds the step size in x, y and z and the size of the grid in x, y and z.
"""

from __future__ import print_function

import numpy as np

class Coordinates(object):
    """Class containing coordinates used in pyrodeo.

    Args:
        x (ndarray): 2D ndarray containing x coordinates
        y (ndarray): 2D ndarray containing y coordinates
        z (ndarray): 2D ndarray containing z coordinates

    Note:
        The validity of the arrays is not checked! To be used in a simulation, they should have the same shape, with x[:,j,k] containing the x coordinates for all j,k, y[i,:,k] containing the y coordinates for all i,k and z[i,j,:] containing the z coordinates for all i,j. In addition, x, y and z should have a constant step size.

    The following public attributes are available:

    Attributes:
        x (ndarray): 2D ndarray containing x coordinates
        y (ndarray): 2D ndarray containing y coordinates
        z (ndarray): 3D ndarray containing z coordinates
        dimensions ([int, int, int]): grid dimensions in the x, y and z direction
        dxyz ([float, float, float]): step size in the x, y and z direction

    """

    def __init__(self, x, y, z):
        self.dimensions = (len(x[:,0,0]), len(y[0,:,0]), len(z[0,0,:]))
        dx = x[1,0,0] - x[0,0,0]
        dy = dx
        dz = dx
        if len(y[0,:,0]) > 1:
            dy = y[0,1,0] - y[0,0,0]
        if len(z[0,0,:]) > 1:
            dz = z[0,0,1] - z[0,0,0]
        self.dxyz = [dx, dy, dz]

        self.x = x
        self.y = y
        self.z = z

    @classmethod
    def from_1d(cls, x, y, z):
        """Initialize from 1D arrays.

        Build coordinates from existing 1D ndarrays (should include ghost cells). Calculate dimensions and step sizes.

        Args:
            x (ndarray): 1D ndarray containing x coordinates
            y (ndarray): 1D ndarray containing y coordinates
            z (ndarray): 1D ndarray containing z coordinates

        Note:
            The validity of arrays x, y and z is not checked. They should have constant step size.

        """
        x, y, z = np.meshgrid(x, y, z, indexing='ij')
        return cls(x, y, z)

    @classmethod
    def from_dims(cls,
                  dimensions=(100, 1, 1),
                  domain=([-0.5, 0.5], [], [])):
        """Initialize from dimensions and domain size.

        Build coordinates given the dimensions of the grid and the size of the domain. Some basic checks are performed to ensure the resulting coordinates are valid.

        Args:
            dimensions (:obj:`(int,int,int)`, optional): Dimensions of the grid
            domain (:obj:`([float,float],[float,float],[float,float])`, optional): Domain boundaries in x, y and z
        """
        # Check if dimensions valid
        if len(dimensions) != 3:
            raise TypeError('Expected dimensions to have three elements')
        if (dimensions[1] != 1 and len(domain[1]) != 2):
            raise TypeError('Invalid y domain')
        if (dimensions[2] != 1 and len(domain[2]) != 2):
            raise TypeError('Invalid z domain')
        if dimensions[0] <= 1:
            raise ValueError('First grid dimension should be larger than 1')

        # Check if domain is valid
        if len(domain) != 3:
            raise TypeError('Expected domain to have three elements')
        if domain[0][1] <= domain[0][0]:
            raise ValueError('Invalid x domain')

        # Step sizes (all x for now, do y and z later if necessary)
        dxy = [(domain[0][1] - domain[0][0])/np.float(dimensions[0]),
               (domain[0][1] - domain[0][0])/np.float(dimensions[0]),
               (domain[0][1] - domain[0][0])/np.float(dimensions[0])]

        # x coordinates, including two ghost cells on either side
        x = np.linspace(domain[0][0] - 1.5*dxy[0],
                        domain[0][1] + 2.5*dxy[0],
                        dimensions[0] + 4,
                        endpoint=False)

        # Add y if necessary
        if dimensions[1] > 1:
            if domain[1][1] <= domain[1][0]:
                raise ValueError('Invalid y domain')
            dxy[1] = (domain[1][1] - domain[1][0])/np.float(dimensions[1])
            y = np.linspace(domain[1][0] - 1.5*dxy[1],
                            domain[1][1] + 2.5*dxy[1],
                            dimensions[1] + 4,
                            endpoint=False)
        else:
            y = [0.0]

        # Add z if necessary
        if dimensions[2] > 1:
            if domain[2][1] <= domain[2][0]:
                raise ValueError('Invalid z domain')
            dxy[2] = (domain[2][1] - domain[2][0])/np.float(dimensions[2])
            z = np.linspace(domain[2][0] - 1.5*dxy[2],
                            domain[2][1] + 2.5*dxy[2],
                            dimensions[2] + 4,
                            endpoint=False)
        else:
            z = [0.0]

        return cls.from_1d(x, y, z)

    def __str__(self):
        """Show information about coordinates."""
        return 'Coordinate object, part of pyrodeo. \n \
    Grid dimensions (nx, ny, nz): {} \n \
    Domain size x (including ghost zones): {}, {} \n \
    Domain size y (including ghost zones): {}, {} \n \
    Domain size z (including ghost zones): {}, {}\n'.format(np.shape(self.x),
                                                            np.min(self.x),
                                                            np.max(self.x),
                                                            np.min(self.y),
                                                            np.max(self.y),
                                                            np.min(self.z),
                                                            np.max(self.z))
