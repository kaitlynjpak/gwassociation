Usage
=====

``gwassociation`` evaluates whether gravitational-wave (GW) sky maps and
electromagnetic (EM) transient candidates are compatible in sky position,
distance, and time.  The package supports two common workflows:

* a GW sky map compared with a point-like EM candidate, and
* a primary GW sky map compared with a secondary sky map.

Installation
------------

Install the package from a checkout with editable mode:

.. code-block:: console

   pip install -e .

Optional extras are available for development, documentation, HEALPix plotting,
and LIGO sky-map integration:

.. code-block:: console

   pip install -e ".[dev,docs,healpy,ligo]"

Point-candidate Python API
--------------------------

The :class:`gwassociation.Association` facade accepts candidate coordinates in
degrees using ``ra`` and ``dec``.  Optional ``z`` and ``z_err`` values activate
the distance-overlap term, and optional ``time``/``gw_time`` values activate the
temporal-overlap term.

.. code-block:: python

   from gwassociation import Association

   assoc = Association(
       "fits_files/S190425z_bayestar.fits.gz,0",
       {
           "name": "AT2024abc",
           "ra": 120.5,
           "dec": -30.0,
           "z": 0.05,
           "z_err": 0.002,
           "time": 1234567890.0,
           "gw_time": 1234567880.0,
       },
   )

   results = assoc.compute_odds(em_model="kilonova")
   print(results["posterior_odds"], results["confidence"])

Candidate coordinate aliases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Candidate dictionaries can also use ``ra_deg``/``dec_deg``,
``ra_rad``/``dec_rad``, or an Astropy ``SkyCoord`` object under the ``coord``
key.  The high-level API normalizes these forms before creating a transient.

.. code-block:: python

   Association(
       "skymap.fits",
       {"ra_deg": 120.5, "dec_deg": -30.0, "time": 1234567890.0},
   )

Secondary-skymap workflow
-------------------------

Provide ``secondary_skymap`` when the second event has a localization map rather
than a point-like sky position.  In this mode, ``compute_odds`` evaluates map
spatial overlap, radial overlap if both maps have distance metadata, and a time
factor based on the two event times.

.. code-block:: python

   assoc = Association(
       "primary.fits",
       secondary_skymap="secondary.fits",
       secondary_event_time=1234567900.0,
   )
   coincidence = assoc.compute_odds(chance_coincidence_rate=1e-4)

Ranking multiple candidates
---------------------------

Use :meth:`gwassociation.association.Association.rank_candidates` when a single
GW map should be compared to several transient dictionaries.

.. code-block:: python

   candidates = [
       {"name": "AT2024abc", "ra": 120.5, "dec": -30.0, "z": 0.05},
       {"name": "AT2024xyz", "ra": 250.1, "dec": 15.4, "z": 0.09},
   ]
   ranked = assoc.rank_candidates(candidates, em_model="kilonova")

Command line
------------

The console script mirrors the point-candidate workflow and writes result files
to the selected output directory.

.. code-block:: console

   gwassociation \
     --gw-file fits_files/S190425z_bayestar.fits.gz,0 \
     --ra 120.5 \
     --dec -30.0 \
     --z 0.05 \
     --time 1234567890 \
     --gw-time 1234567880 \
     --model kilonova \
     --out results

Diagnostic plots
----------------

The high-level API can plot the GW sky map with the EM candidate overlaid when
``healpy`` is installed:

.. code-block:: python

   assoc.plot_skymap("candidate_skymap.png")

Additional helpers in :mod:`gwassociation.plotting.distributions` visualize
approximate distance and temporal distributions for sanity checks.

Interpreting outputs
--------------------

``compute_odds`` returns a dictionary containing the individual overlap terms
(``I_omega``, ``I_dl``, ``I_t``), their Bayes-factor combination, posterior odds,
log-odds, and a confidence-like odds-to-probability conversion.  The absolute
calibration depends on the caller's priors and chance-coincidence rate; inspect
individual terms before comparing results across surveys or counterpart models.
