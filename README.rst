Pyrodeo
========================

.. image:: https://badge.fury.io/py/pyrodeo.svg
    :target: https://badge.fury.io/py/pyrodeo
.. image:: https://readthedocs.org/projects/pip/badge/?version=stable&style=flat
    :target: http://pyrodeo.readthedocs.org
.. image:: http://img.shields.io/badge/license-GPL-green.svg?style=flat
    :target: https://github.com/SijmeJan/pyrodeo/blob/master/LICENSE
.. image:: https://img.shields.io/pypi/pyversions/pyrodeo.svg?style=flat-square
    :target: https://pypi.python.org/pypi/pyrodeo
.. image:: https://img.shields.io/pypi/implementation/pyrodeo.svg?style=flat-square
    :target: https://pypi.python.org/pypi/pyrodeo
.. image:: http://img.shields.io/travis/SijmeJan/pyrodeo/master.svg?style=flat
    :target: https://travis-ci.org/SijmeJan/pyrodeo/
.. image:: https://coveralls.io/repos/SijmeJan/pyrodeo/badge.svg?branch=master&service=github
    :target: https://coveralls.io/github/SijmeJan/pyrodeo?branch=master

Pyrodeo is a Python implementation of the isothermal hydrodynamic
solver RODEO  (ROe solver for Disc-Embedded Objects). Its main purpose
is to perform numerical simulations of astrophysical (gas) discs.

Features
-----------------------------

* Two- and three-dimensional inviscid isothermal hydrodynamics using a
  Riemann solver.
* Second order in space and time in regions of smooth flow.
* Different geometries: Cartesian, Shearing Sheet, cylindrical
  coordinates and spherical coordinates.
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

Changelog
=========

Version 0.0.8
--------------

* Small fixes in tests and documentation

 Version 0.0.7
--------------

* 3D isothermal
* Spherical coordinates
* Logarithmic radial coordinate
* Various small fixes

Version 0.0.6
--------------

* Make `source_func` a source term integrator rather than just returning the extra source terms.
* Move tests inside package.

Version 0.0.5
--------------

* Initial release
