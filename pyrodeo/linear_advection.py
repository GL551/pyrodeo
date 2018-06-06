# -*- coding: utf-8 -*-
"""Defines a linear advection solver in a periodic domain.
"""

from __future__ import print_function

import numpy as np
from pyrodeo.claw_solver import ClawSolver

class LinearAdvection(ClawSolver):
    """Linear advection solver in periodic domain.

    Args:
        advection_velocity (ndarray): Advection velocity. Must have the same shape as :class:`.State` members (density, velocity) to be advected and must be uniform in the first dimension (x).
        limiter_parameter (float): Parameter setting the limiter function. Should be >= 1 (minmod limiter, most diffusive) and <= 2 (superbee limiter, least diffusive).

    Note:
        The linear advection equation is solved on a periodic x domain. The advection velocity can not depend on x. If advecting over y, the state will have to be transposed first.

    The following public attributes and methods are available:

    Attributes:
        advection_velocity (ndarray): Advection velocity. Must have the same shape as :class:`.State` members (density, velocity) to be advected and must be uniform in the first dimension (x).
        sb (float): Parameter setting the limiter function. Should be >= 1 (minmod limiter, most diffusive) and <= 2 (superbee limiter, least diffusive).

    """
    def __init__(self, advection_velocity, limiter_parameter):
        if (limiter_parameter < 1.0 or limiter_parameter > 2.0):
            raise ValueError('Expecting limiter_parameter >= 1 and <= 2 to be TVD')

        self.advection_velocity = advection_velocity
        self.sb = limiter_parameter

    def step(self, dt, dx, state):
        """Perform one time step dt.

        Evolve the linear advection equation over a time step dt, updating the state. The domain can be multidimensional, but evolution is with respect to x. There is no restriction on the magnitude of dt.

        Args:
            dt (float): Time step to take.
            dx (float): Step size in x.
            state (State): :class:`.State` containing density and velocity

        """
        # Work with momenta
        state.velx = state.velx*state.dens
        state.vely = state.vely*state.dens

        # Length of state in x and y
        nx = len(state.dens[:,0])
        ny = len(state.dens[0,:])

        # Advection velocity
        u = self.advection_velocity[2:nx-2,:]

        # Distance to shift solution in x
        shiftx = dt*u

        # Number of integer cells to shift
        Nshift = np.around(shiftx/dx)
        # Remaining u*dt/dx
        udtdx = shiftx/dx - Nshift

        # Views without ghost cells
        dens = state.dens[2:nx-2,:]
        momx = state.velx[2:nx-2,:]
        momy = state.vely[2:nx-2,:]

        # Limited slopes for 2nd order correction
        slopedens = self.limiter(dens - np.roll(dens, 1, axis=0),
                                 np.roll(dens, -1, axis=0) - dens, self.sb)
        slopemomx = self.limiter(momx - np.roll(momx, 1, axis=0),
                                 np.roll(momx, -1, axis=0) - momx, self.sb)
        slopemomy = self.limiter(momy - np.roll(momy, 1, axis=0),
                                 np.roll(momy, -1, axis=0) - momy, self.sb)

        # Upwind state
        s = np.sign(udtdx)
        densu = -0.5*(s - 1.0)*np.roll(dens, -1, axis=0) + \
            0.5*(s + 1.0)*np.roll(dens, 1, axis=0)
        momxu = -0.5*(s - 1.0)*np.roll(momx, -1, axis=0) + \
            0.5*(s + 1.0)*np.roll(momx, 1, axis=0)
        momyu = -0.5*(s - 1.0)*np.roll(momy, -1, axis=0) + \
            0.5*(s + 1.0)*np.roll(momy, 1, axis=0)

        # Upwind slopes
        slopedensu = -0.5*(s - 1.0)*np.roll(slopedens, -1, axis=0) + \
            0.5*(s + 1.0)*np.roll(slopedens, 1, axis=0)
        slopemomxu = -0.5*(s - 1.0)*np.roll(slopemomx, -1, axis=0) + \
            0.5*(s + 1.0)*np.roll(slopemomx, 1, axis=0)
        slopemomyu = -0.5*(s - 1.0)*np.roll(slopemomy, -1, axis=0) + \
            0.5*(s + 1.0)*np.roll(slopemomy, 1, axis=0)

        # Update state
        a = np.abs(udtdx);
        dens += (-a*(dens - densu) - \
                 0.5*udtdx*(1.0 - a)*(slopedens - slopedensu))
        momx += (-a*(momx - momxu) - \
                 0.5*udtdx*(1.0 - a)*(slopemomx - slopemomxu))
        momy += (-a*(momy - momyu) - \
                 0.5*udtdx*(1.0 - a)*(slopemomy - slopemomyu))

        # Shift integer number of cells
        for i in range(ny):
            n = int(Nshift[0,i])
            dens[:,i] = np.roll(dens[:,i], n)
            momx[:,i] = np.roll(momx[:,i], n)
            momy[:,i] = np.roll(momy[:,i], n)

        # Switch back to velocities
        state.velx /= state.dens
        state.vely /= state.dens
