"""Data ingestion helpers for GW-EM association analyses."""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
from typing import Any

from astropy.io import fits
from astropy.table import Table


def _coerce(value: str | None) -> Any:
    if value in (None, "", "None", "null"):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def ingest_transient_list(filename: str | Path) -> list[dict[str, Any]]:
    """Load transient candidates from CSV, JSON, FITS, or whitespace text."""
    path = Path(filename)
    suffix = path.suffix.lower()

    if suffix == ".json":
        return json.loads(path.read_text())

    if suffix == ".csv":
        with path.open(newline="") as handle:
            return [{key: _coerce(value) for key, value in row.items()} for row in csv.DictReader(handle)]

    if suffix in {".txt", ".dat"}:
        lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
        if not lines:
            return []
        header = lines[0].split()
        rows = []
        for line in lines[1:]:
            values = line.split()
            rows.append({key: _coerce(value) for key, value in zip(header, values)})
        return rows

    if suffix in {".fits", ".fit"}:
        table = Table.read(path)
        return [dict(zip(table.colnames, row)) for row in table]

    raise ValueError(f"Unsupported transient-list format: {path}")


def load_gw_posteriors(filename: str | Path) -> dict[str, Any]:
    """Load posterior samples from JSON, HDF5, ASCII, or FITS table files."""
    path = Path(filename)
    suffixes = {suffix.lower() for suffix in path.suffixes}

    if ".json" in suffixes:
        return json.loads(path.read_text())

    if suffixes & {".h5", ".hdf5"}:
        if importlib.util.find_spec("h5py") is None:  # pragma: no cover
            raise ImportError("h5py is required to read HDF5 posterior files.")
        import h5py

        posteriors: dict[str, Any] = {}
        with h5py.File(path, "r") as handle:
            group = handle.get("posterior_samples", handle)
            for key, value in group.items():
                if isinstance(value, h5py.Dataset):
                    posteriors[key] = value[()]
        return posteriors

    if suffixes & {".txt", ".dat"}:
        rows = ingest_transient_list(path)
        if not rows:
            return {}
        return {key: [row.get(key) for row in rows] for key in rows[0]}

    if suffixes & {".fits", ".fit"}:
        table = Table.read(path)
        return {col: table[col].data.tolist() for col in table.colnames}

    raise ValueError(f"Cannot read posterior file: {path}")


def ingest_gracedb_data(superevent_id: str) -> dict[str, Any]:
    """Return GraceDB URL metadata without performing network access."""
    base_url = "https://gracedb.ligo.org/superevents"
    return {
        "superevent_id": superevent_id,
        "skymap_url": f"{base_url}/{superevent_id}/files/bayestar.fits",
        "classification": None,
        "event_time": None,
        "far": None,
        "instruments": None,
        "note": "Metadata stub only; download the file separately for analysis.",
    }


def format_for_analysis(transient_dict: dict[str, Any]) -> dict[str, Any]:
    """Normalize a raw transient dictionary for :class:`gwassociation.Association`."""
    required = ["ra", "dec"]
    optional = ["z", "z_err", "time", "time_err", "magnitude", "filter_band", "name"]
    formatted: dict[str, Any] = {}
    for key in required:
        if key not in transient_dict:
            raise ValueError(f"Missing required field: {key}")
        formatted[key] = float(transient_dict[key])
    for key in optional:
        value = transient_dict.get(key)
        if key in {"z", "z_err", "time", "time_err", "magnitude"} and value is not None:
            formatted[key] = float(value)
        else:
            formatted[key] = value
    return formatted


def load_skymap_metadata(skymap_file: str | Path) -> dict[str, Any]:
    """Load common metadata fields from a FITS sky map header."""
    with fits.open(skymap_file) as hdul:
        header = hdul[1].header if len(hdul) > 1 else hdul[0].header
    return {
        "gps_time": header.get("GPS_TIME"),
        "distance": header.get("DISTMEAN"),
        "distance_std": header.get("DISTSTD"),
        "instruments": header.get("INSTRUME"),
        "creator": header.get("CREATOR", "unknown"),
        "origin": header.get("ORIGIN", "unknown"),
    }


def batch_ingest_transients(file_list: list[str | Path]) -> list[dict[str, Any]]:
    """Load and concatenate transient candidates from multiple files."""
    all_transients: list[dict[str, Any]] = []
    for filename in file_list:
        all_transients.extend(ingest_transient_list(filename))
    return all_transients


def save_results(results: dict[str, Any], output_file: str | Path, format: str = "json") -> None:
    """Save analysis results as JSON or CSV."""
    path = Path(output_file)
    path.parent.mkdir(parents=True, exist_ok=True)
    if format == "json" or path.suffix.lower() == ".json":
        path.write_text(json.dumps(results, indent=2, default=str))
        return
    if format == "csv" or path.suffix.lower() == ".csv":
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=results.keys())
            writer.writeheader()
            writer.writerow(results)
        return
    raise ValueError(f"Unsupported output format: {format}")


__all__ = [
    "ingest_transient_list",
    "load_gw_posteriors",
    "ingest_gracedb_data",
    "format_for_analysis",
    "load_skymap_metadata",
    "batch_ingest_transients",
    "save_results",
]
