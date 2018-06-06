# -*- coding: utf-8 -*-
"""Coordinate class used in pyrodeo.

The Coordinate class hold the x and y coordinates as 2D ndarrays. In addition, it holds the step size in x and y and the size of the grid in x and y.
"""

from __future__ import print_function

import numpy as np

class Coordinates(object):
    """Class containing coordinates used in pyrodeo.

    Args:
        x (ndarray): 2D ndarray containing x coordinates
        y (ndarray): 2D ndarray containing y coordinates

    Note:
        The validity of the arrays is not checked! To be used in a simulation, they should have the same shape, with x[:,i] containing the x coordinates for all i and y[j,:] containing the y coordinates for all j. In addition, both x and y should have a constant step size.

    The following public attributes are available:

    Attributes:
        x (ndarray): 2D ndarray containing x coordinates
        y (ndarray): 2D ndarray containing y coordinates
        dimensions ([int, int]): grid dimensions in x and y direction
        dxy ([float, float]): step size in x and y direction

    """

    def __init__(self, x, y):
        self.dimensions = (len(x[:,0]), len(y[0,:]))
        dx = x[1,0] - x[0,0]
        dy = dx
        if len(y[0,:]) > 1:
            dy = y[0,1] - y[0,0]
        self.dxy = [dx, dy]

        self.x = x
        self.y = y

    @classmethod
    def from_1d(cls, x, y):
        """Initialize from 1D arrays.

        Build coordinates from existing 1D ndarrays (should include ghost cells). Calculate dimensions and step sizes.

        Args:
            x (ndarray): 1D ndarray containing x coordinates
            y (ndarray): 1D ndarray containing y coordinates

        Note:
            The validity of arrays x and y is not checked. They should have constant step size.

        """
        x, y = np.meshgrid(x, y, indexing='ij')
        return cls(x, y)

    @classmethod
    def from_dims(cls,
                  dimensions=(100, 1),
                  domain=([-0.5, 0.5], [])):
        """Initialize from dimensions and domain size.

        Build coordinates given the dimensions of the grid and the size of the domain. Some basic checks are performed to ensure the resulting coordinates are valid.

        Args:
            dimensions (:obj:`(int,int)`, optional): Dimensions of the grid
            domain (:obj:`([float,float],[float,float])`, optional): Domain boundaries in x and y
        """
        # Check if dimensions valid
        if len(dimensions) != 2:
            raise TypeError('Expected dimensions to have two elements')
        if (dimensions[1] != 1 and len(domain[1]) != 2):
            raise TypeError('Invalid y domain')
        if dimensions[0] <= 1:
            raise ValueError('First grid dimension should be larger than 1')

        # Check if domain is valid
        if len(domain) != 2:
            raise TypeError('Expected domain to have two elements')
        if domain[0][1] <= domain[0][0]:
            raise ValueError('Invalid x domain')

        # Step size x and y
        dxy = [(domain[0][1] - domain[0][0])/np.float(dimensions[0]),
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

        return cls.from_1d(x, y)

    def __str__(self):
        """Show information about coordinates."""
        return 'Coordinate object, part of pyrodeo. \n \
    Grid dimensions (nx, ny): {} \n \
    Domain size x (including ghost zones): {}, {} \n \
    Domain size y (including ghost zones): {}, {}\n'.format(np.shape(self.x),
                                                            np.min(self.x),
                                                            np.max(self.x),
                                                            np.min(self.y),
                                                            np.max(self.y))
