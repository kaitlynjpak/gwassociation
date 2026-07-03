"""GW--GW sky-map overlap screening for gravitational-lensing candidates.

This subpackage screens pairs of gravitational-wave sky maps for unusually high
spatial overlap -- a first-pass filter for strongly lensed event pairs (two
detections of a single source) and for prioritising multimessenger follow-up.

The pipeline mirrors the published poster *"Screening for Gravitationally
Lensed Gravitational-Wave Events via Skymap Overlap"* (Pak & Magana Hernandez):

1. :mod:`~gwassociation.screening.gracedb` -- query GraceDB superevents.
2. :mod:`~gwassociation.screening.download` -- fetch the best sky map per event.
3. :mod:`~gwassociation.screening.overlap` -- MOC-aware overlap integral.
4. :mod:`~gwassociation.screening.pipeline` -- all-pairs (parallel) computation.
5. :mod:`~gwassociation.screening.statistics` -- empirical p-values and ranking.
6. :mod:`~gwassociation.screening.plots` -- distribution and joint sky-map plots.

Unlike :class:`gwassociation.association.Association` (which pairs one GW event
with one EM transient), the screening pipeline compares GW sky maps against each
other. It reuses the same Singer et al. (2014) overlap integral for the spatial
term but operates directly on multi-order (MOC/UNIQ) sky maps as served by
GraceDB.
"""

from .overlap import read_skymap, skymap_overlap_integral

__all__ = ["read_skymap", "skymap_overlap_integral"]
