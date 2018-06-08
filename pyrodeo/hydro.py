# -*- coding: utf-8 -*-
"""Defines the Hydro class performing hydrodynamic updates of the state.
"""

from __future__ import print_function

import numpy as np
from pyrodeo.state import State
from pyrodeo.linear_advection import LinearAdvection
from pyrodeo.roe import Roe

class Hydro:
    """Construct class for hydrodynamic updates.

    Construct from existing instances of :class:`.Param` and :class:`.Coordinates`. It constructs a :class:`.LinearAdvection` instance which will deal with orbital advection, and a :class:`.Roe` instance dealing with hydrodynamics of residual velocities.

    Args:
        param (:class:`.Param`): Valid Param object, containing simulation parameters.
        coords (:class:`.Coordinates`): Valid Coordinates object, containing x and y coordinates. Used to calculate orbital advection velocity.

    The following attributes and methods are available:

    Attributes:
        orbital_advection (:class:`.LinearAdvection`): Instance of :class:`.LinearAdvection` class used to do orbital advection.
        roe (:class:`.Roe`): Instance of :class:`.Roe` class dealing with hydrodynamics of residual velocities.

    """
    def __init__(self, param, coords):
        v_adv = 0.0*coords.x
        if param.geometry == 'sheet':
            v_adv = (-1.5*coords.x).transpose()
        if param.geometry == 'cyl':
            v_adv = (np.power(coords.x, -1.5) - param.Omega).transpose()
        self.orbital_advection = LinearAdvection(v_adv, param.limiter_param)
        self.roe = Roe(param.limiter_param, param.min_dens)

    def calc_time_step(self, geometry, coords, state):
        """Calculate time step obeying the CFL condition.

        Args:
            geometry (str): 'cart', 'sheet' or 'cyl'.
            coords (:class:`.Coordinates`): Valid :class:`.Coordinates` object, containing x and y coordinates.
            state (:class:`.State`): Valid :class:`.State` object, containing density and velocity.

        Returns:
            float: Maximum time step obeying the CFL condition.

        """
        dtx = coords.dxy[0]/np.max(np.abs(state.velx) + state.soundspeed)
        dty = dtx
        if len(state.dens[0,:]) > 1:
            cs = state.soundspeed
            if geometry == 'cyl':
                cs = state.soundspeed/coords.x
            dty = coords.dxy[1]/np.max(np.abs(state.vely) + cs)

        return np.min([dtx, dty])

    def shear_periodic_boundaries(self, t, coords, state):
        """Set shear-periodic boundary conditions.

        In the shearing sheet geometry, the x direction can be quasi-periodic, i.e. periodic but modified for the shear. Imagine neighbouring sheets shearing past the center sheet.

        Args:
            t (float): Current simulation time. Used to calculate over what distance neighbouring sheets have shifted since t = 0.
            coords (:class:`.Coordinates`): Valid :class:`.Coordinates` object, containing x and y coordinates.
            state (:class:`.State`): Valid :class:`.State` object, containing density and velocity.
        """
        nx = np.shape(state.dens)[0]

        # Create State made just of ghost zones
        dens = np.concatenate((state.dens[nx-4:nx-2,:], state.dens[2:4,:]))
        velx = np.concatenate((state.velx[nx-4:nx-2,:], state.velx[2:4,:]))
        vely = np.concatenate((state.vely[nx-4:nx-2,:], state.vely[2:4,:]))
        temp_state = State(dens, velx, vely, vely)

        # Advection velocity is speed of next box
        v_adv = 0.0*dens + 1.5*(coords.x[nx-2,0] + coords.x[nx-3,0])
        v_adv[2:4,:] = -v_adv[2:4,:]

        # Linear advection since t = 0 when solution was periodic
        temp_state.transpose()
        v_adv = v_adv.transpose()

        la = LinearAdvection(v_adv, self.orbital_advection.sb)
        la.step(t, coords.dxy[1], temp_state)

        temp_state.transpose()

        # Store in ghost zones
        state.dens[:2,:] = temp_state.dens[:2,:]
        state.velx[:2,:] = temp_state.velx[:2,:]
        state.vely[:2,:] = temp_state.vely[:2,:]
        state.dens[nx-2:,:] = temp_state.dens[2:,:]
        state.velx[nx-2:,:] = temp_state.velx[2:,:]
        state.vely[nx-2:,:] = temp_state.vely[2:,:]

    def preprocess(self, coords, param, state, direction):
        """Modify state to quasi-cartesian form and calculate geometric source terms.

        Isothermal hydrodynamics allows for a generic form of the equations in all geometries, subject only to modifications in the geometric source terms.

        Args:
            coords (:class:`.Coordinates`): Valid :class:`.Coordinates` object, containing x and y coordinates.
            param (:class:`.Param`): Valid :class:`.Param` object, containing simulation parameters.
            state (:class:`.State`): Valid :class:`.State` object, containing density and velocity.
            direction (int): 0 (integrating x) or 1 (integrating y).

        Returns:
            ndarray: Geometric source term of the same shape as state.dens.

        """
        source = 0.0*state.dens

        if param.geometry == 'sheet':
            if direction == 0:
                state.vely += 0.5*param.Omega*coords.x
                source = 2.0*state.dens*param.Omega*(state.vely -
                                                     0.5*param.Omega*coords.x)
        if param.geometry == 'cyl':
            if direction == 0:
                state.dens *= coords.x
                state.vely = coords.x*coords.x*(state.vely +
                                                np.power(coords.x, -1.5))
                source = \
                    (state.vely*state.vely/coords.x -
                     1.0)*state.dens/(coords.x*coords.x) + \
                state.soundspeed*state.soundspeed*state.dens/coords.x
            if direction == 1:
                state.soundspeed /= coords.x

        if direction == 1:
            source = source.transpose()
            state.transpose()
            state.swap_velocities()

        return source

    def postprocess(self, coords, param, state, direction):
        """Inverse of :meth:`.preprocess`.

        Reverse modifications by :meth:`.preprocess`.

        Args:
            coords (:class:`.Coordinates`): Valid :class:`.Coordinates` object, containing x and y coordinates.
            param (:class:`.Param`): Valid :class:`.Param` object, containing simulation parameters.
            state (:class:`.State`): Valid :class:`.State` object, containing density and velocity.
            direction (int): 0 (integrating x) or 1 (integrating y).

        """
        if direction == 1:
            state.transpose()
            state.swap_velocities()

        if param.geometry == 'sheet':
            if direction == 0:
                state.vely -= 0.5*param.Omega*coords.x

        if param.geometry == 'cyl':
            if direction == 0:
                state.dens /= coords.x
                state.vely = state.vely/(coords.x*coords.x) - \
                np.power(coords.x, -1.5)
            if direction == 1:
                state.soundspeed *= coords.x

        if np.min(state.dens) < 0.0:
            i = np.unravel_index(np.argmin(state.dens, axis=None),
                                 state.dens.shape)
            print("Negative density encountered after integrating" +
                  " direction {} at x = {}, y = {}".format(direction,
                                                           coords.x[i],
                                                           coords.y[i]))

    def evolve(self, t, t_max, coords, param, state,
               source_func, source_param):
        """Evolve state from t to tmax.

        Args:
            t (float): Current simulation time.
            t_max (float): Simulation time to reach before stopping.
            coords (:class:`.Coordinates`): Valid :class:`.Coordinates` object, containing x and y coordinates.
            param (:class:`.Param`): Valid :class:`.Param` object, containing simulation parameters.
            state (:class:`.State`): Valid :class:`.State` object, containing density and velocity.
            source_func (callable): Function returning any extra source terms (non-geometric). It should accept the following arguments: t, coords, state, source_param and return density source, x velocity source, y velocity source.
            source_param (array-like): Extra parameters for source_func.

        Returns:
            (tuple): Tuple consisting of:

                t (float): new simulation time (= t_max if no problems encountered).

                :class:`.State`: Updated :class:`.State`.

        """
        # Evolve until t = t_max
        while (t < t_max):
            # Calculate time step
            dt = param.courant*self.calc_time_step(param.geometry,
                                                   coords, state)
            # End exactly on t_max
            if (t + dt > t_max):
                dt = t_max - t

            # Copy state in case something goes wrong
            old_state = State.copy(state)

            # Continue until everything OK
            success = False
            while success is False:
                print("Current time: {}, time step: {}".format(t, dt))

                # Shear periodic x boundaries if necessary
                if (param.geometry == 'sheet' and
                    param.boundaries[0] == 'shear periodic'):
                    self.shear_periodic_boundaries(t, coords, state)

                # Dimensional split: do all dimensions independently
                for dim in (0,1):
                    if np.shape(state.dens)[dim] > 1:
                        source = self.preprocess(coords, param, state, dim)

                        # Hydrodynamic update
                        self.roe.step(dt, coords.dxy[dim], state,
                                      source, param.boundaries[dim])

                        # Orbital advection
                        if dim == 1:
                            self.orbital_advection.step(dt,
                                                        coords.dxy[dim],
                                                        state)

                        self.postprocess(coords, param, state, dim)

                # Check if density positive
                success = True
                mindens = np.min(state.dens)

                # If not, restore state and try with smaller time step
                if mindens < 0.0:
                    state = State.copy(old_state)
                    dt = 0.5*dt
                    success = False

            # Integrate extra source terms (ignoring ghost zones!)
            if source_func is not None:
                sdens, svelx, svely = \
                  source_func(t, coords, state, source_param)
                state.dens += dt*sdens*state.no_ghost
                state.velx += dt*svelx*state.no_ghost
                state.vely += dt*svely*state.no_ghost

            # Update simulation time
            t = t + dt

        return t, state
