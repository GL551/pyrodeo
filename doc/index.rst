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

Equations solved
-------------------------------------

The current version supports inviscid isothermal hydrodynamics in two spatial
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
    v_x) + \frac{\partial}{\partial y}(\rho v_y)=0

Momentum in x-direction:

.. math::

    \frac{\partial}{\partial t}(\rho v_x) + \frac{\partial}{\partial x}(\rho
    v_x^2 + c^2\rho) + \frac{\partial}{\partial y}(\rho v_x v_y)=0

Momentum in y-direction:

.. math::

    \frac{\partial}{\partial t}(\rho v_y) + \frac{\partial}{\partial x}(\rho
    v_x v_y) + \frac{\partial}{\partial y}(\rho v_y^2 + c^2\rho)=0

Shearing Sheet
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Shearing Sheet is essentially a Cartesian model of a small patch
in an astrophysical disc. This patch is rotating at the local
Keplerian velocity, which means that Coriolis and centrifugal-type forces need to be
included on the right-hand side of the equations. On the other hand,
the patch is assumed to be small enough so that a local Cartesian
frame can be used in stead of cylindrical coordinates. We therefore
still have the continuity equation:

.. math::

    \frac{\partial\rho}{\partial t} + \frac{\partial}{\partial x}(\rho
    v_x) + \frac{\partial}{\partial y}(\rho v_y)=0

The x-momentum equation now includes source terms on the right-hand side:

.. math::

    \frac{\partial}{\partial t}(\rho v_x) + \frac{\partial}{\partial x}(\rho
    v_x^2 + c^2\rho) + \frac{\partial}{\partial y}(\rho v_x
    v_y)=2\Omega\rho v_y + 3\rho\Omega^2 x

Same for the momentum equation in y-direction:

.. math::

    \frac{\partial}{\partial t}(\rho v_y) + \frac{\partial}{\partial x}(\rho
    v_x v_y) + \frac{\partial}{\partial y}(\rho v_y^2 +
    c^2\rho)=-2\Omega\rho v_x

Cylindrical coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a full disc cylindrical coordinates are preferred. This time we
have geometrical source terms and gravity from the central object to
worry about. The continuity equation now reads:

.. math::

    \frac{\partial}{\partial t}(r\rho) + \frac{\partial}{\partial r}(r\rho
    v_r) + \frac{\partial}{\partial \varphi}(r\rho v_\varphi)=0

Note that :math:`v_\varphi` is an angular velocity. The radial
momentum equation now includes source terms representing centrifugal
and gravitational forces, in addition to a geometrical pressure source
term:

.. math::

    \frac{\partial}{\partial t}(r\rho v_r) + \frac{\partial}{\partial r}(r\rho
    v_r^2 + c^2r\rho) + \frac{\partial}{\partial \varphi}(\rho v_r
    v_\varphi)= r^2\rho v_\varphi^2 - \rho\frac{GM_*}{r} + c^2\rho

In the :math:`\varphi` direction we get a Coriolis source term:

.. math::

    \frac{\partial}{\partial t}(r\rho v_\varphi) + \frac{\partial}{\partial r}(r\rho
    v_r v_\varphi) + \frac{\partial}{\partial \varphi}(r\rho v_\varphi^2 +
    c^2\rho/r)=-2\rho v_r v_\varphi

Numerical method
-------------------------------

Dimensional splitting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pyrodeo uses dimensional splitting to integrate the equations. For the
x direction (therefore neglecting y-derivatives), we can cast the
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

For the Cartesian setup, we simply have

.. math::

   \bar\rho =\rho, \bar v_x = v_x, \bar v_y = v_y, \bar c = c, \bar x
   = x, \bar y = y, S_x = 0

For the shearing sheet, we need

.. math::

   \bar\rho =\rho, \bar v_x = v_x, \bar v_y = v_y - \Omega x/2, \bar c = c, \bar x
   = x, \bar y = y, S_x = 2\Omega \rho (v_y +3\Omega x/2)

Finally, in cylindrical coordinates we need

.. math::

   \bar\rho =r\rho, \bar v_x = v_r, \bar v_y = r^2v_\varphi, \bar c = c, \bar x
   = r, \bar y = \varphi, S_x = r^2\rho v_\varphi^2 -
   \rho\frac{GM_*}{r} + c^2\rho

For the y-integration (neglecting x-derivatives) we can cast the
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

For both the Cartesian setup and the shearing sheet, we simply have

.. math::

   \bar\rho =\rho, \bar v_x = v_x, \bar v_y = v_y, \bar c = c, \bar x
   = x, \bar y = y, S_y = 0

In cylindrical coordinates we need

.. math::

   \bar\rho =\rho, \bar v_x = v_r, \bar v_y = v_\varphi, \bar c = c/r, \bar x
   = r, \bar y = \varphi, S_y = 0

Note that the resulting equations are very symmetric in x and y: if we
swap x and y in the y integration the equations have exactly the same
form as for the x integration. Therefore, if we prepare all quantities
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



Examples
-------------------------------

A simple example in Cartesian geometry is an isothermal shock tube:

.. literalinclude:: ../examples/test_cartesian.py

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

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
