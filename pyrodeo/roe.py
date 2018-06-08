# -*- coding: utf-8 -*-
"""Defines the Roe class defining the Roe solver.
"""

from __future__ import print_function

import numpy as np
from pyrodeo.claw_solver import ClawSolver

class Roe(ClawSolver):
    """Construct class for the Roe solver.

    Constructor sets two basic attributes (sb and min_dens), after which a time step can be taken through the :meth:`.step` method.

    Args:
        flux_limiter (float): Flux limiter parameter. Should be >= 1 (minmod, most diffusive limiter) and <= 2 (superbee, least diffusive limiter).
        min_dens (float): Minimum density when to switch to HLL to preserve positivity.

    The following attributes and methods are available:

    Attributes:
        sb (float): Flux limiter parameter. Should be >= 1 (minmod, most diffusive limiter) and <= 2 (superbee, least diffusive limiter).
        min_dens (float): Minimum density when to switch to HLL to preserve positivity.
    """
    def __init__(self, flux_limiter, min_dens):
        self.sb = flux_limiter
        self.min_dens = min_dens

    def limit_flux(self, dens, dens_left, f1dens, f2dens, dtdx):
        """Limit second order flux to preserve positivity

        Calculate the maximum contribution of the second order flux in order for the density to remain positive.

        Args:
            dens (ndarray): Current density.
            dens_left (ndarray): Density in the cell to the left
            f1dens (ndarray): First order mass flux.
            f2dens (ndarray): Second order mass flux.
            dtdx (float): Time step / space step

        Returns:
            ndarray: Array with values >= 0 and <= 1 specifying the maximum contribution of second order flux for the density to remain positive.

        """
        thetap = 0.0*dens + 1.0
        rhop = 0.5*(dens_left - 2.0*dtdx*f2dens)
        epsd = np.minimum(0.0*dens + self.min_dens,
                          np.minimum(dens_left, dens))
        sel = np.where(rhop < epsd)
        rho1 = 0.5*(dens_left[sel] - 2.0*dtdx*f1dens[sel])
        if len(rho1) >= 1:
            if np.min(rho1) < 0.0:
                dens[sel] = -1.0*dens[sel]
        thetap[sel] = np.maximum(0.0*rho1,
                                 (rho1 - epsd[sel])/(rho1 - rhop[sel]))
        if np.max(thetap) > 1.0:
            dens = -1.0*dens
        thetam = 0.0*dens + 1.0
        rhom = 0.5*(dens + 2.0*dtdx*f2dens)
        sel = np.where(rhom < epsd)
        rho1 = 0.5*(dens[sel] + 2.0*dtdx*f1dens[sel])
        if len(rho1) >= 1:
            if np.min(rho1) < 0.0:
                dens[sel] = -1.0*dens[sel]
        thetam[sel] = np.maximum(0.0*rho1,
                                 (rho1 - epsd[sel])/(rho1 - rhom[sel]))
        if np.max(thetam) > 1.0:
            dens = -1.0*dens

        return np.minimum(thetap, thetam)

    def step(self, dt, dx, state, source, bc):
        """Update state for a single time step.

        Args:
            dt (float): Time step.
            dx (float): Space step.
            state (:class:`.State`): Current :class:`.State`, will be updated.
            source (ndarray): Geometric source terms, must have same shape as state.dens.
            bc (str): Boundary condition: 'periodic' or 'reflect' (other boundary conditions are dealt with elsewhere).

        """
        #nx = len(state.dens[:,0])
        #ny = len(state.dens[0,:])

        # Set periodic boundaries
        if bc == 'periodic':
            #state.dens[:2,:] = state.dens[nx-4:nx-2,:]
            #state.velx[:2,:] = state.velx[nx-4:nx-2,:]
            #state.vely[:2,:] = state.vely[nx-4:nx-2,:]
            state.dens[:2,:] = state.dens[-4:-2,:]
            state.velx[:2,:] = state.velx[-4:-2,:]
            state.vely[:2,:] = state.vely[-4:-2,:]

            #state.dens[nx-2:,:] = state.dens[2:4,:]
            #state.velx[nx-2:,:] = state.velx[2:4,:]
            #state.vely[nx-2:,:] = state.vely[2:4,:]
            state.dens[-2:,:] = state.dens[2:4,:]
            state.velx[-2:,:] = state.velx[2:4,:]
            state.vely[-2:,:] = state.vely[2:4,:]

        dens_left = np.roll(state.dens, 1, axis=0)
        velx_left = np.roll(state.velx, 1, axis=0)
        vely_left = np.roll(state.vely, 1, axis=0)

        cs_left = np.roll(state.soundspeed, 1, axis = 0)
        # Average sound speed 0.5*((i) + (i-1))
        cs_left = np.roll(state.soundspeed, 1, axis = 0)
        cs = 0.5*(state.soundspeed + cs_left)

        dtdx = dt/dx

        # Parameter vector
        z0 = np.sqrt(state.dens)
        z1 = z0*state.velx
        z2 = z0*state.vely

        # Average parameter vector: 0.5*((i) + (i-1))
        z0 = 0.5*(z0 + np.roll(z0, 1, axis=0))
        z1 = 0.5*(z1 + np.roll(z1, 1, axis=0))
        z2 = 0.5*(z2 + np.roll(z2, 1, axis=0))

        # Roe averages
        u = z1/z0
        v = z2/z0

        # HLL signal speeds
        bm = np.minimum(0.0*u, np.minimum(u - cs,
                                          velx_left - cs_left))
        bp = np.maximum(0.0*u, np.maximum(u + cs,
                                          state.velx + state.soundspeed))

        # Cell-centred fluxes
        fdens = state.dens*state.velx
        fmomx = fdens*state.velx + \
          state.soundspeed*state.soundspeed*state.dens
        fmomy = fdens*state.vely

        # Left and right fluxes (note stationary extrapolation)
        fldens = np.roll(fdens, 1, axis=0)
        flmomx = np.roll(fmomx + 0.5*dx*source, 1, axis=0)
        flmomy = np.roll(fmomy, 1, axis=0)
        frdens = fdens
        frmomx = fmomx - 0.5*dx*source
        frmomy = fmomy

        # Flux difference: (i) - (i-1)
        dfdens =  frdens - fldens
        dfmomx =  frmomx - flmomx
        dfmomy =  frmomy - flmomy

        # State difference: (i) - (i-1)
        ddens = state.dens - dens_left
        dmomx = state.dens*state.velx - dens_left*velx_left
        dmomy = state.dens*state.vely - dens_left*vely_left

        # HLL flux
        f1dens = (bp*fldens - bm*frdens + bp*bm*ddens)/(bp - bm);
        f1momx = (bp*flmomx - bm*frmomx + bp*bm*dmomx)/(bp - bm);
        f1momy = (bp*flmomy - bm*frmomy + bp*bm*dmomy)/(bp - bm);

        # Projection coefficients
        b1 =-0.5*(dfmomx - (u + cs)*dfdens)/cs
        b2 = 0.5*(dfmomx - (u - cs)*dfdens)/cs
        b3 = dfmomy - v*dfdens;

        b1 = b1/(u - cs + 1.0e-10)
        b2 = b2/(u + cs + 1.0e-10)
        b3 = b3/(u + 1.0e-10)

        # Upwind projection coefficients
        b1u = -0.5*(np.sign(u - cs) - 1.0)*np.roll(b1, -1, axis=0) + \
          0.5*(np.sign(u - cs) + 1.0)*np.roll(b1, 1, axis=0)
        b2u = -0.5*(np.sign(u + cs) - 1.0)*np.roll(b2, -1, axis=0) + \
          0.5*(np.sign(u + cs) + 1.0)*np.roll(b2, 1, axis=0)
        b3u = -0.5*(np.sign(u) - 1.0)*np.roll(b3, -1, axis=0) + \
          0.5*(np.sign(u) + 1.0)*np.roll(b3, 1, axis=0)

        # Flux limiter
        fl1 = self.limiter(b1, b1u, self.sb)
        fl2 = self.limiter(b2, b2u, self.sb)
        fl3 = self.limiter(b3, b3u, self.sb)

        # Second order projection coefficients
        b1 = -np.abs(u - cs)*b1 + fl1*(np.abs(u - cs) - dtdx*(u - cs)*(u - cs))
        b2 = -np.abs(u + cs)*b2 + fl2*(np.abs(u + cs) - dtdx*(u + cs)*(u + cs))
        b3 = -np.abs(u)*b3 + fl3*(np.abs(u) - dtdx*u*u)

        # Interface flux (i - 1/2)
        f2dens = 0.5*(fldens + frdens + b1 + b2)
        f2momx = 0.5*(flmomx + frmomx + (u - cs)*b1 + (u + cs)*b2)
        f2momy = 0.5*(flmomy + frmomy + v*(b1 + b2) + b3)

        # Limit flux to maintain positivity
        theta = self.limit_flux(state.dens, dens_left, f1dens, f2dens, dtdx)

        # Limited fluxes
        fdens = theta*f2dens - (theta - 1.0)*f1dens
        fmomx = theta*f2momx - (theta - 1.0)*f1momx
        fmomy = theta*f2momy - (theta - 1.0)*f1momy

        if bc == 'reflect':
            fdens[2,:] = 0.0
            fmomx[2,:] = state.soundspeed[2,:]*state.soundspeed[2,:]* \
              state.dens[2,:] - 0.5*dx*source[2,:]
            fmomy[2,:] = 0.0
            fdens[-2,:] = 0.0
            fmomx[-2,:] = state.soundspeed[-3,:]* \
              state.soundspeed[-3,:]*state.dens[-3,:] + \
              0.5*dx*source[-3,:]
            fmomy[-2,:] = 0.0

        # Interface flux difference: (i+1) - (i)
        dfdens = np.roll(fdens, -1, axis=0) - fdens
        dfmomx = np.roll(fmomx, -1, axis=0) - fmomx
        dfmomy = np.roll(fmomy, -1, axis=0) - fmomy

        # Change in density and momenta (ignoring ghost cells)
        ddens = -dtdx*dfdens*state.no_ghost
        dmomx = (dt*source - dtdx*dfmomx)*state.no_ghost
        dmomy = -dtdx*dfmomy*state.no_ghost

        state.dens += ddens
        state.velx += (dmomx - ddens*state.velx)/state.dens
        state.vely += (dmomy - ddens*state.vely)/state.dens

        ## # Ghost cells are unaffected
        ## if ny > 1:
        ##     ddens[:,0] = 0.0
        ##     dmomx[:,0] = 0.0
        ##     dmomy[:,0] = 0.0

        ##     ddens[:,1] = 0.0
        ##     dmomx[:,1] = 0.0
        ##     dmomy[:,1] = 0.0

        ##     ddens[:,ny-1] = 0.0
        ##     dmomx[:,ny-1] = 0.0
        ##     dmomy[:,ny-1] = 0.0

        ##     ddens[:,ny-2] = 0.0
        ##     dmomx[:,ny-2] = 0.0
        ##     dmomy[:,ny-2] = 0.0

        ## # View without ghost cells x
        ## dens = state.dens[2:nx-2,:]
        ## velx = state.velx[2:nx-2,:]
        ## vely = state.vely[2:nx-2,:]

        ## dens += ddens[2:nx-2,:]
        ## velx += (dmomx[2:nx-2,:] - ddens[2:nx-2,:]*velx)/dens
        ## vely += (dmomy[2:nx-2,:] - ddens[2:nx-2,:]*vely)/dens
