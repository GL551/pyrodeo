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
        coords (:class:`.Coordinates`): Valid Coordinates object, containing x, y and z coordinates. Used to calculate orbital advection velocity.

    The following attributes and methods are available:

    Attributes:
        orbital_advection (:class:`.LinearAdvection`): Instance of :class:`.LinearAdvection` class used to do orbital advection.
        roe (:class:`.Roe`): Instance of :class:`.Roe` class dealing with hydrodynamics of residual velocities.

    """
    def __init__(self, param, coords):
        v_adv = 0.0*coords.x
        if param.geometry == 'sheet':
            v_adv = (-1.5*coords.x).transpose((1,0,2))
        if param.geometry == 'cyl':
            v_adv = (np.power(coords.x, -1.5) -
                     param.frame_rotation).transpose((1,0,2))
        if param.geometry == 'sph':
            v_adv = (np.power(coords.x, -1.5)/np.sin(coords.z) -
                     param.frame_rotation).transpose((1,0,2))
        self.orbital_advection = LinearAdvection(v_adv, param.limiter_param)
        self.roe = Roe(param.limiter_param, param.min_dens)

    def calc_time_step(self, geometry, coords, state, log_radial=False):
        """Calculate time step obeying the CFL condition.

        Args:
            geometry (str): 'cart', 'sheet' or 'cyl'.
            coords (:class:`.Coordinates`): Valid :class:`.Coordinates` object, containing x, y and z coordinates.
            state (:class:`.State`): Valid :class:`.State` object, containing density and velocity.
            log_radial (:obj:`bool`, optional): Flag indicating whether a logarithmic radial coordinate is used

        Returns:
            float: Maximum time step obeying the CFL condition.

        """
        cs = state.soundspeed

        dtx = 1.0e10
        dty = dtx
        dtz = dtx
        if len(state.dens[:,0,0]) > 1:
            abs_speed = np.abs(state.velx) + cs
            if log_radial is True:
                abs_speed /= coords.x
            dtx = coords.dxyz[0]/np.max(abs_speed)

        if len(state.dens[0,:,0]) > 1:
            abs_speed = np.abs(state.vely) + cs
            if geometry == 'cyl':
                abs_speed /= coords.x
            if geometry == 'sph':
                abs_speed /= (coords.x*np.sin(coords.z))
            dty = coords.dxyz[1]/np.max(abs_speed)

        if len(state.dens[0,0,:]) > 1:
            abs_speed = np.abs(state.velz) + cs
            if geometry == 'sph':
                abs_speed /= coords.x
            dtz = coords.dxyz[2]/np.max(abs_speed)

        return np.min([dtx, dty, dtz])

    def shear_periodic_boundaries(self, t, coords, state):
        """Set shear-periodic boundary conditions.

        In the shearing sheet geometry, the x direction can be quasi-periodic, i.e. periodic but modified for the shear. Imagine neighbouring sheets shearing past the center sheet.

        Args:
            t (float): Current simulation time. Used to calculate over what distance neighbouring sheets have shifted since t = 0.
            coords (:class:`.Coordinates`): Valid :class:`.Coordinates` object, containing x and y coordinates.
            state (:class:`.State`): Valid :class:`.State` object, containing density and velocity.
        """
        # Single cell in x; nothing to do
        if len(state.dens[:,0,0]) <= 1:
            return

        # Create State made just of ghost zones
        dens = np.concatenate((state.dens[-4:-2,:,:], state.dens[2:4,:,:]))
        velx = np.concatenate((state.velx[-4:-2,:,:], state.velx[2:4,:,:]))
        vely = np.concatenate((state.vely[-4:-2,:,:], state.vely[2:4,:,:]))
        velz = np.concatenate((state.velz[-4:-2,:,:], state.velz[2:4,:,:]))
        temp_state = State(dens, velx, vely, velz, velz)

        # Advection velocity is speed of next box
        v_adv = 0.0*dens + 1.5*(coords.x[-2:-1,0,0] + coords.x[-3:-2,0,0])
        v_adv[2:4,:] = -v_adv[2:4,:,:]

        # Linear advection since t = 0 when solution was periodic
        temp_state.transpose((1,0,2))
        v_adv = v_adv.transpose((1,0,2))

        la = LinearAdvection(v_adv, self.orbital_advection.sb)
        la.step(t, coords.dxyz[1], temp_state)

        temp_state.transpose((1,0,2))

        # Store in ghost zones
        state.dens[:2,:,:]  = temp_state.dens[:2,:,:]
        state.velx[:2,:,:]  = temp_state.velx[:2,:,:]
        state.vely[:2,:,:]  = temp_state.vely[:2,:,:]
        state.velz[:2,:,:]  = temp_state.velz[:2,:,:]

        state.dens[-2:,:,:] = temp_state.dens[2:,:,:]
        state.velx[-2:,:,:] = temp_state.velx[2:,:,:]
        state.vely[-2:,:,:] = temp_state.vely[2:,:,:]
        state.velz[-2:,:,:] = temp_state.velz[2:,:,:]

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
                state.vely += 0.5*param.frame_rotation*coords.x
                source = 2.0*state.dens*param.frame_rotation*\
                  (state.vely - 0.5*param.frame_rotation*coords.x)
            if direction == 2:
                source = -state.dens*param.frame_rotation*\
                  param.frame_rotation*coords.z

        if param.geometry == 'cyl':
            if direction == 0:
                state.dens *= coords.x
                state.vely = coords.x*state.vely + np.sqrt(coords.x)
                dpot = coords.x*np.power(coords.x*coords.x +
                                         coords.z*coords.z, -1.5)
                source = \
                    state.dens*state.vely*state.vely/ \
                    (coords.x*coords.x*coords.x) - \
                    state.dens*dpot + \
                state.soundspeed*state.soundspeed*state.dens/coords.x
                if param.log_radial is True:
                    state.velx = state.velx/coords.x
                    state.dens *= coords.x
                    state.soundspeed /= coords.x
                    source = source - state.dens*state.velx*state.velx - \
                      state.soundspeed*state.soundspeed*state.dens
            if direction == 1:
                state.vely /= coords.x
                state.soundspeed /= coords.x
            if direction == 2:
                dpot = coords.z*np.power(coords.x*coords.x +
                                         coords.z*coords.z, -1.5)
                source = -state.dens*dpot

        if param.geometry == 'sph':
            if direction == 0:
                c = state.soundspeed

                state.dens *= coords.x*coords.x
                state.vely = state.vely*coords.x + np.sqrt(coords.x)
                state.velz *= coords.x
                source = state.dens*(((state.vely*state.vely +
                                       state.velz*state.velz)/coords.x -
                                       1.0)/(coords.x*coords.x) +
                                       2.0*c*c/coords.x)
                if param.log_radial is True:
                    state.velx = state.velx/coords.x
                    state.dens *= coords.x
                    state.soundspeed /= coords.x
                    source = source - state.dens*state.velx*state.velx - \
                      state.soundspeed*state.soundspeed*state.dens
            if direction == 1:
                g = 1.0/(coords.x*np.sin(coords.z))
                state.vely *= g
                state.soundspeed *= g

            if direction == 2:
                sint = np.sin(coords.z)

                state.dens *= sint
                state.vely = sint*(state.vely/coords.x +
                                   np.power(coords.x, -1.5))
                state.velz /= coords.x
                state.soundspeed /= coords.x

                cott = np.cos(coords.z)/sint
                source = state.dens*cott*(state.vely*state.vely/(sint*sint) +
                                          state.soundspeed*state.soundspeed)

        if direction == 1:
            source = source.transpose((1,0,2))
            state.transpose((1,0,2))
            state.swap_velocities(1)

        if direction == 2:
            source = source.transpose((2,1,0))
            state.transpose((2,1,0))
            state.swap_velocities(2)

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
            state.transpose((1,0,2))
            state.swap_velocities(1)

        if direction == 2:
            state.transpose((2,1,0))
            state.swap_velocities(2)

        if param.geometry == 'sheet':
            if direction == 0:
                state.vely -= 0.5*param.frame_rotation*coords.x

        if param.geometry == 'cyl':
            if direction == 0:
                state.dens /= coords.x
                state.vely = state.vely/coords.x - np.power(coords.x, -0.5)
                if param.log_radial is True:
                    state.velx *= coords.x
                    state.dens /= coords.x
                    state.soundspeed *= coords.x

            if direction == 1:
                state.vely *= coords.x
                state.soundspeed *= coords.x

        if param.geometry == 'sph':
            if direction == 0:
                state.dens /= (coords.x*coords.x)
                state.vely = (state.vely - np.sqrt(coords.x))/coords.x
                state.velz /= coords.x
                if param.log_radial is True:
                    state.velx *= coords.x
                    state.dens /= coords.x
                    state.soundspeed *= coords.x

            if direction == 1:
                g = coords.x*np.sin(coords.z)
                state.vely *= g
                state.soundspeed *= g

            if direction == 2:
                sint = np.sin(coords.z)
                state.dens /= sint
                state.vely = state.vely*coords.x/sint - \
                  np.power(coords.x, -0.5)
                state.velz *= coords.x
                state.soundspeed *= coords.x

    def evolve(self, t, t_max, coords, param, state,
               source_func, source_param):
        """Evolve state from t to tmax.

        Args:
            t (float): Current simulation time.
            t_max (float): Simulation time to reach before stopping.
            coords (:class:`.Coordinates`): Valid :class:`.Coordinates` object, containing x and y coordinates.
            param (:class:`.Param`): Valid :class:`.Param` object, containing simulation parameters.
            state (:class:`.State`): Valid :class:`.State` object, containing density and velocity.
            source_func (callable): Function integrating any extra source terms (non-geometric). It should accept the following arguments: t, dt, coords, state, source_param.
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
                                                   coords, state,
                                                   param.log_radial)
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
                for dim in (0,1,2):
                    if np.shape(state.dens)[dim] > 1:
                        source = self.preprocess(coords, param, state, dim)

                        # Hydrodynamic update
                        self.roe.step(dt, coords.dxyz[dim], state,
                                      source, param.boundaries[dim])

                        # Orbital advection
                        if dim == 1:
                            self.orbital_advection.step(dt,
                                                        coords.dxyz[dim],
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
                source_func(t, dt, coords, state, source_param)

            # Update simulation time
            t = t + dt

        return t, state
