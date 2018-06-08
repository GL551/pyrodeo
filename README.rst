Pyrodeo
========================

.. image:: http://img.shields.io/badge/license-GPL-green.svg?style=flat
    :target: https://github.com/SijmeJan/pyrodeo/blob/master/LICENSE

Pyrodeo is a Python implementation of the isothermal hydrodynamic
solver RODEO  (ROe solver for Disc-Embedded Objects). Its main purpose
is to perform numerical simulations of astrophysical (gas) discs.

Features
-----------------------------

* Two-dimensional inviscid isothermal hydrodynamics using a Riemann
  solver.
* Second order in space and time in regions of smooth flow.
* Different geometries: Cartesian, Shearing Sheet and cylindrical
  coordinates.
* HDF5 output

Quick start
-----------------------------

Pyrodeo can be installed from the command line simply by entering::

  pip install pyrodeo

A simple simulation can be created and run by entering:

.. code:: python

          import pyrodeo
          sim = pyrodeo.Simulation.from_geom('cart')
          sim.evolve([0.25])

Since the standard initial conditions consist of constant density and
pressure and zero velocity, no visible evolution takes place. For more
interesting examples, see the documentation.

Documentation
-------------
The full documentation can be found at

http://pyrodeo.readthedocs.org
