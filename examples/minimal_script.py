"""Minimal gwassociation API example.

Run from the repository root with a FITS sky map available locally, e.g.::

    python examples/minimal_script.py fits_files/S190425z_bayestar.fits.gz,0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running directly from a source checkout without installing first.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gwassociation import Association


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a minimal GW-EM association calculation.")
    parser.add_argument("gw_file", help="Path to a 2D or 3D GW sky-map FITS file.")
    parser.add_argument("--ra", type=float, default=120.5, help="Transient RA in degrees.")
    parser.add_argument("--dec", type=float, default=-30.2, help="Transient Dec in degrees.")
    parser.add_argument("--z", type=float, default=0.05, help="Transient redshift.")
    parser.add_argument("--time", type=float, default=1234567890.0, help="Transient event time.")
    parser.add_argument("--gw-time", type=float, default=1234567880.0, help="GW trigger time.")
    args = parser.parse_args()

    assoc = Association(
        args.gw_file,
        {
            "name": "example_candidate",
            "ra": args.ra,
            "dec": args.dec,
            "z": args.z,
            "time": args.time,
            "gw_time": args.gw_time,
        },
    )
    results = assoc.compute_odds()
    print(f"P(associated): {results['confidence']:.3%}")
    print(f"Posterior odds: {results['posterior_odds']:.3e}")


if __name__ == "__main__":
    main()
