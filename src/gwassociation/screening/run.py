"""End-to-end driver for the overlap screening pipeline.

:func:`run_screening` chains the stages -- query, download, all-pairs overlap,
statistics, ranking, and plotting -- and writes the data products to an output
directory. It is the single entry point used by the ``gwassociation screen``
command, and is equally usable from Python.
"""

from __future__ import annotations

import json
import os
from dataclasses import asdict

from . import download as download_mod
from . import gracedb, pipeline, plots, statistics


def _default_progress():
    """Return tqdm if available, else an identity progress wrapper."""
    try:
        from tqdm import tqdm

        return tqdm
    except Exception:  # pragma: no cover - tqdm is part of the screen extra
        def _identity(iterable, **_kwargs):
            return iterable

        return _identity


def run_screening(
    out_dir: str,
    skymap_dir: str | None = None,
    time_window: str = gracedb.DEFAULT_TIME_WINDOW,
    far_threshold: float = gracedb.DEFAULT_FAR_THRESHOLD,
    bbh_threshold: float | None = gracedb.DEFAULT_BBH_THRESHOLD,
    workers: int = 1,
    top_n: int = 30,
    min_days_apart: int = 0,
    n_plot_pairs: int = 1,
    reuse_overlaps: str | None = None,
    log=print,
) -> dict:
    """Run the screening pipeline and write data products to ``out_dir``.

    Parameters
    ----------
    out_dir : str
        Output directory for ``overlaps.json``, ``top_pairs.csv``, ``summary.json``
        and the figures.
    skymap_dir : str, optional
        Directory for the sky-map download cache. Defaults to ``out_dir/skymaps``.
    time_window, far_threshold, bbh_threshold :
        GraceDB query cuts (see :func:`gracedb.query_superevents`).
    workers : int
        Worker processes for the all-pairs overlap computation.
    top_n : int
        Number of ranked pairs to write to CSV.
    min_days_apart : int
        Drop pairs whose events are fewer than this many days apart (use 1 to
        exclude same-day duplicate detections).
    n_plot_pairs : int
        How many of the top-ranked pairs to render joint sky-map plots for.
    reuse_overlaps : str, optional
        Path to a previously saved overlaps JSON. When given, the query,
        download, and overlap stages are skipped and these values are used --
        handy for re-ranking or re-plotting without recomputing.
    log : callable
        Sink for progress messages (``print`` by default; pass a no-op to mute).

    Returns
    -------
    dict
        The summary that is also written to ``summary.json``.
    """
    os.makedirs(out_dir, exist_ok=True)
    if skymap_dir is None:
        skymap_dir = os.path.join(out_dir, "skymaps")
    progress = _default_progress()

    overlaps_path = os.path.join(out_dir, "overlaps.json")

    if reuse_overlaps:
        log(f"Loading precomputed overlaps from {reuse_overlaps}")
        overlaps = pipeline.load_overlaps(reuse_overlaps)
        query_info = None
        download_info = None
    else:
        # --- Stage 1: query GraceDB ---
        log(f"Querying GraceDB: created {time_window!r}, far < {far_threshold}")
        result = gracedb.query_superevents(
            time_window=time_window,
            far_threshold=far_threshold,
            bbh_threshold=bbh_threshold,
            progress=progress,
        )
        log(f"  {len(result.event_ids)} events kept "
            f"({result.n_kept_unclassified} lacked classification)")
        query_info = asdict(result)

        if len(result.event_ids) < 2:
            log(
                f"WARNING: the GraceDB query matched {len(result.event_ids)} event(s); "
                "at least 2 are needed to form a pair. Outputs will be empty. "
                f"Check that --time-window ({time_window!r}) covers a period with "
                "public superevents and that --far-threshold is not too strict."
            )

        # --- Stage 2: download sky maps ---
        log(f"Downloading sky maps to {skymap_dir}")
        dl = download_mod.download_skymaps(
            result.event_ids, skymap_dir, progress=progress
        )
        log(f"  downloaded {len(dl.downloaded)}, failed {len(dl.failed)}")
        download_info = {"n_downloaded": len(dl.downloaded), "n_failed": len(dl.failed)}

        # --- Stage 3: all-pairs overlap ---
        n_expected = pipeline.expected_pair_count(skymap_dir)
        log(f"Computing overlaps for {n_expected:,} pairs using {workers} worker(s)")
        overlaps = pipeline.compute_overlaps(skymap_dir, workers=workers, progress=progress)
        pipeline.save_overlaps(overlaps, overlaps_path)
        log(f"  saved overlaps to {overlaps_path}")

    # --- Stage 4: statistics ---
    stats = statistics.describe(overlaps)
    pvalues = statistics.empirical_pvalues(overlaps)
    corrections = statistics.multiple_comparison_corrections(pvalues)
    log(f"Valid pairs: {stats.get('n', 0):,}; "
        f"surviving Bonferroni: {corrections.n_bonferroni}, FDR: {corrections.n_fdr}")

    # --- Stage 5: ranking + CSV ---
    ranked = statistics.rank_pairs(
        overlaps, pvalues=pvalues, top_n=top_n, min_days_apart=min_days_apart
    )
    csv_path = os.path.join(out_dir, "top_pairs.csv")
    statistics.write_top_pairs_csv(ranked, csv_path)
    log(f"Wrote top {len(ranked)} pairs to {csv_path}")

    # --- Stage 6: plots ---
    figures = {}
    try:
        hist_path = os.path.join(out_dir, "overlap_histogram.png")
        plots.plot_overlap_distribution(overlaps, hist_path)
        figures["histogram"] = hist_path

        surv_path = os.path.join(out_dir, "survival_function.png")
        plots.plot_survival_function(overlaps, surv_path)
        figures["survival"] = surv_path
    except Exception as exc:  # pragma: no cover - plotting is best-effort
        log(f"Warning: distribution plots failed: {exc}")

    # Pairs to plot: the top n_plot_pairs from the ranking, plus the highest-
    # overlap pair whose events are at least one day apart. The overall top pair
    # is often a same-day duplicate detection of a single trigger; the top
    # >=1-day pair is the genuine time-separated (lensing) candidate, so it is
    # always plotted even when it is not in the leading n_plot_pairs.
    plot_rows = list(ranked[:n_plot_pairs])
    distinct = statistics.rank_pairs(overlaps, pvalues=pvalues, top_n=1, min_days_apart=1)
    top_distinct = distinct[0] if distinct else None
    if top_distinct and not any(
        r["event1"] == top_distinct["event1"] and r["event2"] == top_distinct["event2"]
        for r in plot_rows
    ):
        log(f"Also plotting top >=1-day-apart pair: "
            f"{top_distinct['event1']} & {top_distinct['event2']} "
            f"(O={top_distinct['overlap']:.2f}, {top_distinct['days_apart']} d)")
        plot_rows.append(top_distinct)

    joint_areas = _plot_pairs(plot_rows, skymap_dir, out_dir, log)
    figures.update(joint_areas["figures"])

    # --- Summary ---
    summary = {
        "parameters": {
            "time_window": time_window,
            "far_threshold": far_threshold,
            "bbh_threshold": bbh_threshold,
            "workers": workers,
            "top_n": top_n,
            "min_days_apart": min_days_apart,
            "skymap_dir": skymap_dir,
            "reuse_overlaps": reuse_overlaps,
        },
        "query": query_info,
        "download": download_info,
        "statistics": stats,
        "corrections": asdict(corrections),
        "top_pairs": ranked,
        "top_distinct_day_pair": top_distinct,
        "joint_areas": joint_areas["areas"],
        "outputs": {
            "overlaps_json": overlaps_path if not reuse_overlaps else reuse_overlaps,
            "top_pairs_csv": csv_path,
            "figures": figures,
        },
    }
    with open(os.path.join(out_dir, "summary.json"), "w") as handle:
        json.dump(summary, handle, indent=2)
    log(f"Wrote summary to {os.path.join(out_dir, 'summary.json')}")
    return summary


def _plot_pairs(rows, skymap_dir, out_dir, log):
    """Render joint sky-map plots for an explicit list of pair rows."""
    figures = {}
    areas = {}
    for row in rows:
        e1, e2 = row["event1"], row["event2"]
        try:
            p1 = download_mod.find_skymap_file(e1, skymap_dir)
            p2 = download_mod.find_skymap_file(e2, skymap_dir)
        except FileNotFoundError:
            log(f"Warning: sky maps for {e1} & {e2} not in {skymap_dir}; skipping plot")
            continue
        save_path = os.path.join(out_dir, f"joint_{e1}_{e2}.png")
        try:
            areas[f"{e1},{e2}"] = plots.plot_joint_skymap(
                p1, p2, save_path,
                event1_id=e1, event2_id=e2, overlap_value=row["overlap"],
            )
            figures[f"joint_{e1}_{e2}"] = save_path
        except Exception as exc:  # pragma: no cover - plotting is best-effort
            log(f"Warning: joint sky-map plot for {e1} & {e2} failed: {exc}")
    return {"figures": figures, "areas": areas}
