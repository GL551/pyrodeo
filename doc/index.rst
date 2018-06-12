.. pyrodeo documentation master file, created by
   sphinx-quickstart on Fri Jun  1 21:12:37 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyrodeo!
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Pyrodeo is a Python implementation of RODEO (ROe solver for
Disc-Embedded Objects), a Roe solver implementation aimed at
hydrodynamic simulations of astrophysical discs.

Installation
-------------------------------------

If Python is not installed, download from `here
<https://www.python.org/downloads/>`_ and install. The latest versions
of Python come with package manager pip included. Then Pyrodeo can be
installed from the command line simply by entering::

  pip install pyrodeo

Quick start
-------------------------------------

Within Python, first import the simulation module::

  >>> import pyrodeo

Create a simulation in Cartesian geometry with standard parameters::

  >>> sim = pyrodeo.Simulation.from_geom('cart')

Run the simulation up to t=0.25::

  >>> sim.evolve([0.25])

Since the standard initial conditions consist of constant density and
pressure and zero velocity, no visible evolution takes place. For more
interesting examples, see :ref:`Examples`.


Equations solved
-------------------------------------

The current version supports inviscid isothermal hydrodynamics in three spatial
dimensions. Isothermal means the pressure :math:`p` is related to the
density :math:`\rho` simply through :math:`p=c^2\rho`, where the sound
speed :math:`c` is either a constant (fully isothermal) or a
prescribed function of position (locally isothermal). Three geometries
are available: Cartesian coordinates, the Shearing Sheet and
Cylindrical coordinates.

Cartesian coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In Cartesian coordinates, we have the continuity equation:

.. math::

    \frac{\partial\rho}{\partial t} + \frac{\partial}{\partial x}(\rho
    v_x) + \frac{\partial}{\partial y}(\rho v_y) + \frac{\partial}{\partial z}(\rho v_z)=0

Momentum in x-direction:

.. math::

    \frac{\partial}{\partial t}(\rho v_x) + \frac{\partial}{\partial x}(\rho
    v_x^2 + c^2\rho) + \frac{\partial}{\partial y}(\rho v_x v_y) + \frac{\partial}{\partial z}(\rho v_x v_z)=0

Momentum in y-direction:

.. math::

    \frac{\partial}{\partial t}(\rho v_y) + \frac{\partial}{\partial x}(\rho
    v_x v_y) + \frac{\partial}{\partial y}(\rho v_y^2 + c^2\rho) + \frac{\partial}{\partial z}(\rho
    v_y v_z)=0

Momentum in z-direction:

.. math::

    \frac{\partial}{\partial t}(\rho v_z) + \frac{\partial}{\partial x}(\rho
    v_x v_z) + \frac{\partial}{\partial y}(\rho v_y v_z)+ \frac{\partial}{\partial z}(\rho v_z^2 + c^2\rho) =0

Shearing Sheet
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Shearing Sheet is essentially a Cartesian model of a small patch
in an astrophysical disc. This patch is rotating at the local
Keplerian velocity, which means that Coriolis and centrifugal-type forces need to be
included on the right-hand side of the equations. On the other hand,
the patch is assumed to be small enough so that a local Cartesian
frame can be used in stead of cylindrical coordinates. Usually the
computational domain is taken to be periodic in y and shear-periodic
in x (periodic but corrected for the shear). We therefore still have
the continuity equation:

.. math::

    \frac{\partial\rho}{\partial t} + \frac{\partial}{\partial x}(\rho
    v_x) + \frac{\partial}{\partial y}(\rho v_y) + \frac{\partial}{\partial z}(\rho v_z)=0

The x-momentum equation now includes source terms on the right-hand side:

.. math::

    \frac{\partial}{\partial t}(\rho v_x) + \frac{\partial}{\partial x}(\rho
    v_x^2 + c^2\rho) + \frac{\partial}{\partial y}(\rho v_x
    v_y) + \frac{\partial}{\partial z}(\rho v_x v_z)=2\Omega\rho v_y + 3\rho\Omega^2 x

Same for the momentum equation in y-direction:

.. math::

    \frac{\partial}{\partial t}(\rho v_y) + \frac{\partial}{\partial x}(\rho
    v_x v_y) + \frac{\partial}{\partial y}(\rho v_y^2 +
    c^2\rho) + \frac{\partial}{\partial z}(\rho
    v_y v_z)=-2\Omega\rho v_x

In the z-direction we get a source term due to the vertical component
of the stellar gravity:

.. math::

    \frac{\partial}{\partial t}(\rho v_z) + \frac{\partial}{\partial x}(\rho
    v_x v_z) + \frac{\partial}{\partial y}(\rho v_y v_z)+
    \frac{\partial}{\partial z}(\rho v_z^2 + c^2\rho) =-\rho \Omega^2 z

.. NOTE::
   In the shearing sheet the sound speed should really be constant (no
   locally isothermal shearing sheet). Together, sound speed and
   angular velocity define a length scale :math:`c/\Omega`, which is a
   measure of the scale height of the disc. Typically one chooses
   :math:`c=\Omega=1`, so that distances are measured in scale heights
   and time in inverse orbital frequency.


Cylindrical coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a full disc in two dimensions cylindrical coordinates
:math:`(R,\varphi,z)` are preferred. This time we have geometrical
source terms and gravity from the central object to worry about. The
continuity equation now reads:

.. math::

    \frac{\partial}{\partial t}(R\rho) + \frac{\partial}{\partial R}(R\rho
    v_R) + \frac{\partial}{\partial \varphi}(\rho v_\varphi) + \frac{\partial}{\partial z}(R\rho v_z)=0

The radial momentum equation now includes source terms representing
centrifugal and gravitational forces, in addition to a geometrical
pressure source term:

.. math::

    \frac{\partial}{\partial t}(R\rho v_R) + \frac{\partial}{\partial R}(R\rho
    v_R^2 + c^2R\rho) + \frac{\partial}{\partial \varphi}(\rho v_R
    v_\varphi)+\frac{\partial}{\partial z}(R\rho v_R v_z)= \rho v_\varphi^2 - R^2\rho\frac{GM_*}{(R^2+z^2)^{3/2}} + c^2\rho

In the :math:`\varphi` direction we get a Coriolis source term:

.. math::

    \frac{\partial}{\partial t}(\rho v_\varphi) + \frac{\partial}{\partial R}(\rho
    v_R v_\varphi) + \frac{1}{R}\frac{\partial}{\partial \varphi}(\rho v_\varphi^2 +
    c^2\rho) +\frac{\partial}{\partial z}(\rho v_\varphi v_z)=-2\rho v_R v_\varphi/R

In the vertical direction we again have the vertical component of the
stellar gravity:

.. math::

    \frac{\partial}{\partial t}(R\rho v_z) + \frac{\partial}{\partial R}(R\rho
    v_R v_z) + \frac{\partial}{\partial \varphi}(\rho v_\varphi v_z)+
    \frac{\partial}{\partial z}(R\rho v_z^2 + c^2R\rho) = -R\rho z\frac{GM_*}{(R^2+z^2)^{3/2}}

.. NOTE::
   The unit of mass is taken to be the mass of the central object. The
   unit of distance is some reference radius. The unit of time is the
   inverse angular velocity at the reference radius. In this system of
   units, the gravitational constant is unity, and one orbit equals
   :math:`2\pi` time units.

Spherical coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a full disc in three dimensions spherical coordinates
:math:`(r,\theta,\varphi)` are often preferred. The continuity and
momentum equations now read:

.. math::

    \frac{\partial}{\partial t}(r^2\sin\theta\rho) +
    \frac{\partial}{\partial r}(r^2\sin\theta\rho v_r) + \frac{\partial}{\partial \varphi}(r\rho v_\varphi) + \frac{\partial}{\partial \theta}(r\sin\theta\rho v_\theta)=0

.. math::

    \frac{\partial}{\partial t} (r^2\sin\theta\rho v_r)
    +\frac{\partial}{\partial r}(r^2\sin\theta\rho(v_r^2+c^2)) +
    \frac{\partial}{\partial\theta}(r\rho v_r v_\theta\sin\theta) + \frac{\partial}{\partial\varphi} (r\rho v_r v_\varphi) = \\
    r^2\sin\theta\rho\frac{v_\theta^2 + v_\varphi^2}{r}
    -r^2\sin\theta\rho\frac{\partial \Phi}{\partial r}+2r\sin\theta c^2\rho

.. math::

   \frac{\partial}{\partial t} (r^2\sin\theta\rho v_\theta) +
   \frac{\partial}{\partial r}(r^2\sin\theta\rho v_rv_\theta) +
   \frac{\partial}{\partial\theta}(r\sin\theta(\rho v_\theta^2 + p)) +
   \frac{\partial}{\partial\varphi}(r\rho v_\varphi v_\theta) = \\
   r\rho v_\varphi^2\cos\theta - r\sin\theta\rho v_rv_\theta + r\cos\theta p

.. math::

   \frac{\partial}{\partial t} (r^2\sin\theta\rho v_\varphi) +
   \frac{\partial}{\partial r}(r^2\sin\theta\rho v_r v_\varphi) +
   \frac{\partial}{\partial\theta}(r\sin\theta\rho v_\theta
   v_\varphi) + \frac{\partial}{\partial\varphi}(r\rho
   (v_\varphi^2+c^2)) = \\
   -r\sin\theta\rho v_\varphi v_r-r\rho v_\varphi
   v_\theta \cos\theta


Extra source terms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pyrodeo solves inviscid isothermal hydrodynamics, and in the shearing
sheet and  cylindrical geometries only gravity from the central object is
considered. Extra physics, as far as it concerns extra source terms,
can be added by a user-defined source integration function. See the :ref:`Examples`
section.

Numerical method
-------------------------------

Dimensional splitting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pyrodeo uses dimensional splitting to integrate the equations.

First direction (x, r)
"""""""""""""""""""""""""""

For the x direction (therefore neglecting y- and z-derivatives), we can cast the
equations into the following form:

.. math::

    \frac{\partial\bar\rho}{\partial t} + \frac{\partial}{\partial
    \bar x}(\bar\rho\bar v_x) = 0

.. math::

    \frac{\partial}{\partial t}(\bar\rho \bar v_x) +
    \frac{\partial}{\partial \bar x}(\bar\rho
    \bar v_x^2 + \bar c^2\bar \rho) = S_x

.. math::

    \frac{\partial}{\partial t}(\bar\rho \bar v_y) +
    \frac{\partial}{\partial \bar x}(\bar\rho\bar v_x \bar v_y) =0

.. math::

    \frac{\partial}{\partial t}(\bar\rho \bar v_z) +
    \frac{\partial}{\partial \bar x}(\bar\rho\bar v_x \bar v_z) =0

For the Cartesian setup, we simply have

.. math::

   \bar\rho =\rho, \bar v_x = v_x, \bar v_y = v_y, \bar v_z = v_z, \bar c = c, \bar x
   = x, \bar y = y, \bar z = z, S_x = 0

For the shearing sheet, we need

.. math::

   \bar\rho =\rho, \bar v_x = v_x, \bar v_y = v_y - \Omega x/2, \bar v_z = v_z, \bar c = c, \bar x
   = x, \bar y = y, \bar z = z , S_x = 2\Omega \rho (v_y +3\Omega x/2)

In cylindrical coordinates we need

.. math::

   \bar\rho =R\rho, \bar v_x = v_R, \bar v_y = Rv_\varphi, \bar v_z
   = v_z, \bar c = c, \bar x
   = R, \bar y = \varphi, \bar z = z , S_x = \rho v_\varphi^2 -
   \frac{R^2\rho GM_*}{(R^2+z^2)^{3/2}} + c^2\rho

Finally, for spherical coordinates we need:

.. math::

   \bar\rho =r^2\sin\theta \rho, \bar v_x =v_r, \bar v_y =rv_\varphi ,
   \bar v_z =rv_\theta, \bar c = c, \bar x = r, \bar y = \varphi, \bar
   z = \theta, S_x = \bar\rho\frac{\bar v_z^2 + \bar
   v_y^2}{r^3}-\bar\rho\frac{\partial \Phi}{\partial r}+\frac{2
   c^2\bar\rho}{r}

Second direction (y, :math:`\varphi`)
"""""""""""""""""""""""""""""""""""""""""""""""

For the y-integration (neglecting x- and z-derivatives) we can cast the
equations in the form:

.. math::

    \frac{\partial\bar\rho}{\partial t} + \frac{\partial}{\partial
    \bar y}(\bar\rho\bar v_y) = 0

.. math::

    \frac{\partial}{\partial t}(\bar\rho \bar v_x) +
    \frac{\partial}{\partial \bar y}(\bar\rho\bar v_x \bar v_y) =0

.. math::

    \frac{\partial}{\partial t}(\bar\rho \bar v_y) +
    \frac{\partial}{\partial \bar y}(\bar\rho
    \bar v_y^2 + \bar c^2\bar \rho) = S_y

.. math::

    \frac{\partial}{\partial t}(\bar\rho \bar v_z) +
    \frac{\partial}{\partial \bar y}(\bar\rho\bar v_y \bar v_z) =0

For both the Cartesian setup and the shearing sheet, we simply have

.. math::

   \bar\rho =\rho, \bar v_x = v_x, \bar v_y = v_y, \bar v_z = v_z, \bar c = c, \bar x
   = x, \bar y = y, \bar z = z, S_y = 0

In cylindrical coordinates we need

.. math::

   \bar\rho =\rho, \bar v_x = v_R, \bar v_y = v_\varphi/R, \bar v_z =
   v_z, \bar c = c/R, \bar x = R, \bar y = \varphi, \bar z = z, S_y = 0

Finally, spherical coordinates:

.. math::

   \bar\rho =rho, \bar v_x =v_r, \bar v_y =v_\varphi/(r\sin\theta),
   \bar v_z =v_\theta, \bar c =c/(r\sin\theta), \bar x =r, \bar y =\varphi, \bar
   z =\theta, S_y =0

Third direction (z, :math:`\theta`)
"""""""""""""""""""""""""""""""""""""""""""""""

Finally, for the z integration we can cast the equations in the form:

.. math::

    \frac{\partial\bar\rho}{\partial t} + \frac{\partial}{\partial
    \bar z}(\bar\rho\bar v_z) = 0

.. math::

    \frac{\partial}{\partial t}(\bar\rho \bar v_x) +
    \frac{\partial}{\partial \bar y}(\bar\rho\bar v_x \bar v_z) =0

.. math::

    \frac{\partial}{\partial t}(\bar\rho \bar v_y) +
    \frac{\partial}{\partial \bar z}(\bar\rho\bar v_y \bar v_z) =0

.. math::

    \frac{\partial}{\partial t}(\bar\rho \bar v_z) +
    \frac{\partial}{\partial \bar z}(\bar\rho
    \bar v_z^2 + \bar c^2\bar \rho) = S_z

For both the Cartesian setup and the shearing sheet, we simply have

.. math::

   \bar\rho =\rho, \bar v_x = v_x, \bar v_y = v_y, \bar v_z = v_z, \bar c = c, \bar x
   = x, \bar y = y, \bar z = z, S_z = -\rho\Omega^2 z

In cylindrical coordinates we need

.. math::

   \bar\rho =\rho, \bar v_x = v_R, \bar v_y = v_\varphi/R, \bar v_z =
   v_z, \bar c = c, \bar x = R, \bar y = \varphi, \bar z = z, S_z =-\rho z\frac{GM_*}{(R^2+z^2)^{3/2}}

Finally, spherical coordinates:

.. math::

   \bar\rho =\sin\theta \rho, \bar v_x =v_r, \bar v_y =\sin\theta v_\varphi,
   \bar v_z =v_\theta/r, \bar c =c/r, \bar x =r, \bar y =\varphi, \bar
   z =\theta, S_z =\rho \cos\theta (v_\varphi^2+c^2)/r^2

Unified approach
"""""""""""""""""""""""""""

Note that the resulting equations are very symmetric in x, y and z: if we
swap x and y in the y integration the equations have exactly the same
form as for the x integration. Similar for the z integration when
swapping z and x. Therefore, if we prepare all quantities
appropriately, we only need a single hydrodynamic solver that is able
to advance the system

.. math::

    \frac{\partial\bar\rho}{\partial t} + \frac{\partial}{\partial
    \bar x}(\bar\rho\bar v_x) = 0

.. math::

    \frac{\partial}{\partial t}(\bar\rho \bar v_x) +
    \frac{\partial}{\partial \bar x}(\bar\rho
    \bar v_x^2 + \bar c^2\bar \rho) = S_x

.. math::

    \frac{\partial}{\partial t}(\bar\rho \bar v_y) +
    \frac{\partial}{\partial \bar x}(\bar\rho\bar v_x \bar v_y) =0

This is what is done in the :class:`.Roe` class. The necessary
preparation is done in the :class:`.Hydro` class.

Orbital advection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the case of the shearing sheet and the cylindrical disc there is a
large :math:`\bar v_y` that severely limits the time step. This limit
can be overcome by splitting the y integration one more by writing
:math:`\bar v_y = \bar u_y + v_0` where :math:`v_0` is independent of y

.. math::

    \frac{\partial\bar\rho}{\partial t} +
    v_0\frac{\partial\bar\rho}{\partial\bar y} + \frac{\partial}{\partial
    \bar y}(\bar\rho\bar u_y) = 0

.. math::

    \frac{\partial}{\partial t}(\bar\rho \bar v_x) +
    v_0\frac{\partial}{\partial\bar y} (\bar\rho \bar v_x)+
    \frac{\partial}{\partial \bar y}(\bar\rho\bar v_x \bar u_y) =0

.. math::

    \frac{\partial}{\partial t}(\bar\rho \bar u_y) +
    v_0\frac{\partial}{\partial\bar y} (\bar\rho \bar u_y)+
    \frac{\partial}{\partial \bar y}(\bar\rho
    \bar u_y^2 + \bar c^2\bar \rho) = 0

The terms involving :math:`v_0` make up a linear advection problem
that can be solved straightforwardly for any time step. This is done
in the :class:`.LinearAdvection` class. The remaining terms are
integrated in the :class:`.Roe` class, but at a much larger time step
because presumably :math:`u_y \ll v_0`.

Algorithm overview
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A single time step in :func:`Hydro.evolve <pyrodeo.hydro.Hydro.evolve>` consists of the following steps:

1. Calculate time step using :func:`Hydro.calc_time_step
   <pyrodeo.hydro.Hydro.calc_time_step>`.

2. Set shear periodic boundary conditions if necessary.

3. Preprocessing step to cast the equations in the same form for all
   geometries and directions, while at the same time calculating the
   source term using :func:`Hydro.preprocess
   <pyrodeo.hydro.Hydro.preprocess>`.

4. Use the Roe solver to advance the hydrodynamic equations using
   :func:`Roe.step <pyrodeo.roe.Roe.step>`.

5. Do orbital advection if necessary using :func:`Hydro.orbital_advection
   <pyrodeo.hydro.Hydro.orbital_advection>`.

6. Do the inverse of step 3, getting all quantities back to their
   original form in :func:`Hydro.postprocess
   <pyrodeo.hydro.Hydro.postprocess>`.

7. Integrate any extra source terms.

Output
-------------------------------

Once the integration routine :func:`Simulation.evolve
<pyrodeo.simulation.Simulation.evolve>` has finished the final state
is available through :func:`simulation.state`. In addition, an output
file `rodeo.h5` is created containing the state at all specified
checkpoints. This is an `HDF5
<https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_ file
created with `h5py <https://www.h5py.org/>`_. It contains the
following groups:

* param: Simulation parameters as specified in the
  :class:`Param <pyrodeo.param.Param>` class.

* coords: Coordinates from the :class:`Coordinates
  <pyrodeo.coords.Coordinates>` class.

* checkpoint#: State at checkpoint, where # stands for an integer.

.. NOTE::
   Both state and coordinate arrays include two ghost zones on each
   side in all directions. This is in order to be able to restore a
   simulation from a checkpoint.

.. NOTE::
    The value stored in `state.vely` is the y-velocity with the
    orbital advection velocity removed! In other words, the
    equilibrium solution in a constant pressure shearing sheet or
    cylindrical disc has vanishing `state.vely`.

An example of reading the file and plotting using `matplotlib <https://matplotlib.org/>`_:

.. literalinclude:: ../examples/example_plot.py

.. _Examples:

Examples
-------------------------------

Shock tube
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A simple example in Cartesian geometry is a one-dimensional isothermal shock tube:

.. literalinclude:: ../examples/example_cartesian.py

The standard grid dimensions are `(100,1)`, which means 100 cells in x
and 1 in y. Try a higher resolution by explicitly specifying the dimensions in
:func:`Simulation.from_geom <pyrodeo.simulation.Simulation.from_geom>`.

Instability in shearing sheet
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A more demanding two-dimensional calculation involves the instability
of a sharp density ridge in the shearing sheet:

.. literalinclude:: ../examples/example_sheet.py

The checkpoints have been chosen close enough together to allow for
the results to be animated:

.. literalinclude:: ../examples/example_animation.py

Disc-planet interaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As the final, most complex example, consider a planet embedded in a
disc in cylindrical coordinates. Since pyrodeo only includes gravity
from the central star, we need to provide an extra source term to
account for the gravitational force due to the planet. In addition, we
define wave-killing zones on the radial edges of the domain to avoid
wave reflection. It will take some time to run this simulation, so
have a cup of tea and come back to see a Jupiter-like planet carve out
a gap in the disc.

.. literalinclude:: ../examples/example_planet.py

Class reference
=========================

.. automodule:: pyrodeo

Coordinates
-------------------------------

.. automodule:: pyrodeo.coords
   :members:

State
-------------------------------

.. automodule:: pyrodeo.state
   :members:

Param
-------------------------------

.. automodule:: pyrodeo.param
   :members:

Conservation law solver
-------------------------------

.. automodule:: pyrodeo.claw_solver
   :members:

Linear advection solver
-------------------------------

.. automodule:: pyrodeo.linear_advection
   :members:

Roe solver
-------------------------------

.. automodule:: pyrodeo.roe
   :members:

Hydro
-------------------------------

.. automodule:: pyrodeo.hydro
   :members:

Simulation
-------------------------------

.. automodule:: pyrodeo.simulation
   :members:
