Usage
=====

.. _installation:

Installation
------------

Install the package from a checkout with editable mode:

.. code-block:: console

   pip install -e .

Python API
----------

.. code-block:: python

   from gwassociation import Association

   assoc = Association("skymap.fits", {
       "ra": 120.5,
       "dec": -30.0,
       "z": 0.05,
       "time": 1234567890.0,
       "gw_time": 1234567880.0,
   })
   results = assoc.compute_odds()

Command line
------------

.. code-block:: console

   gwassociation --gw-file skymap.fits --ra 120.5 --dec -30.0 --time 1234567890
