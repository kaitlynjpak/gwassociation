"""Command-line interface for ``gwassociation``.

The console script is a command group with two sub-commands:

* ``gwassociation odds`` -- GW--EM association: take a primary sky map and either
  a point-like transient or a secondary sky map and compute posterior odds.
* ``gwassociation screen`` -- GW--GW lensing screening: query GraceDB, download
  sky maps, compute all-pairs overlap integrals, and rank high-overlap pairs.
"""

import json
import pathlib
import click
import numpy as np


@click.group()
@click.version_option(package_name="gwassociation", message="%(version)s")
def main():
    """gwassociation: GW--EM association and GW--GW overlap screening."""


@main.command()
@click.option("--gw-file", type=click.Path(exists=True, dir_okay=False), required=True,
              help="Path to primary GW skymap file (FITS).")
@click.option("--secondary-skymap", type=click.Path(exists=True, dir_okay=False), default=None,
              help="Optional secondary skymap (e.g., EM localization) for coincidence analysis.")
@click.option("--ra", type=float, default=None, help="Transient RA [deg] (required without --secondary-skymap).")
@click.option("--dec", type=float, default=None, help="Transient Dec [deg] (required without --secondary-skymap).")
@click.option("--z", type=float, default=None, help="Transient redshift.")
@click.option("--z-err", type=float, default=None, help="Redshift uncertainty.")
@click.option("--time", "ttime", type=float, default=None,
              help="Transient time (GPS or MJD). Required without --secondary-skymap.")
@click.option("--secondary-time", type=float, default=None,
              help="Event time for the secondary skymap (GPS).")
@click.option("--gw-time", type=float, default=None, help="GW event time (GPS).")
@click.option("--model", type=click.Choice(['kilonova', 'grb', 'afterglow']), 
              default='kilonova', help="EM counterpart model.")
@click.option("--out", "outdir", type=click.Path(file_okay=False), default="out",
              help="Output directory for results.")
@click.option("--verbose", is_flag=True, help="Verbose output.")
def odds(gw_file, secondary_skymap, ra, dec, z, z_err, ttime,
         secondary_time, gw_time, model, outdir, verbose):
    """Run GW-EM association analysis"""
    
    # Only import here to avoid issues if package not fully installed
    from gwassociation import Association
    from gwassociation.plots import plot_association_summary
    
    out = pathlib.Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    
    # Validate inputs
    if secondary_skymap is None:
        missing = [name for name, value in (("RA", ra), ("Dec", dec), ("time", ttime)) if value is None]
        if missing:
            raise click.UsageError(
                f"{', '.join(missing)} required unless --secondary-skymap is provided."
            )
    
    # Set GW time to slightly before transient time if not provided
    if gw_time is None:
        if ttime is not None:
            gw_time = ttime - 86400  # Default to 1 day before transient
        else:
            gw_time = 0.0
    
    transient_payload = None
    if ra is not None and dec is not None and ttime is not None:
        transient_payload = {
            'ra': ra,
            'dec': dec,
            'z': z,
            'z_err': z_err,
            'time': ttime,
            'gw_time': gw_time
        }
    
    assoc = Association(
        gw_file,
        transient_payload,
        secondary_skymap=secondary_skymap,
        secondary_event_time=secondary_time
    )
    
    # Compute odds with proper model
    results = assoc.compute_odds(em_model=model)
    
    if verbose:
        print("\n=== GW-EM Association Analysis Results ===")
        print(f"Primary GW File: {gw_file}")
        if secondary_skymap:
            print(f"Secondary skymap: {secondary_skymap}")
            if secondary_time is not None:
                print(f"Secondary event time: {secondary_time}")
        if transient_payload:
            print(f"Transient: RA={ra:.3f}°, Dec={dec:.3f}°")
            if z is not None:
                print(f"Redshift: z={z:.4f} ± {z_err:.4f}" if z_err else f"Redshift: z={z:.4f}")
            print(f"Time: {ttime:.2f} (transient), {gw_time:.2f} (GW)")
            print(f"EM Model: {model}")
        print(f"\nOverlap Integrals:")
        print(f"  Spatial (I_Ω):  {results['I_omega']:.3e}")
        print(f"  Distance (I_DL): {results['I_dl']:.3e}")
        print(f"  Temporal (I_t):  {results['I_t']:.3e}")
        print(f"\nStatistics:")
        print(f"  Bayes Factor:    {results['bayes_factor']:.3e}")
        print(f"  Posterior Odds:  {results['posterior_odds']:.3e}")
        print(f"  Log₁₀ Odds:      {results['log_posterior_odds']:.2f}")
        print(f"  P(Associated):   {results['confidence']:.1%}")
        print(f"\nDecision: {'ASSOCIATED' if results['associated'] else 'NOT ASSOCIATED'}")
    else:
        print(f"P(Associated) = {results['confidence']:.1%}")
        print(f"Decision: {'ASSOCIATED' if results['associated'] else 'NOT ASSOCIATED'}")
    
    # Generate plots
    try:
        if transient_payload:
            fig_path = out / "skymap.png"
            assoc.plot_skymap(str(fig_path))
            if verbose:
                print(f"\nSaved skymap: {fig_path}")
        else:
            fig_path = None
        
        # Additional plots if we have the plotting module
        try:
            summary_path = out / "association_summary.png"
            plot_association_summary(results, str(summary_path))
            if verbose:
                print(f"Saved summary: {summary_path}")
        except Exception as e:
            if verbose:
                print(f"Warning: Could not generate summary plot: {e}")
                
    except Exception as e:
        if verbose:
            print(f"Warning: Could not generate plots: {e}")
    
    # Save results to JSON
    with open(out / "results.json", "w") as f:
        # Convert numpy types to native Python types for JSON serialization
        json_results = {}
        for k, v in results.items():
            if isinstance(v, (np.float32, np.float64)):
                json_results[k] = float(v)
            elif isinstance(v, (np.int32, np.int64)):
                json_results[k] = int(v)
            elif isinstance(v, np.ndarray):
                json_results[k] = v.tolist()
            elif v == np.inf:
                json_results[k] = "inf"
            elif v == -np.inf:
                json_results[k] = "-inf"
            elif isinstance(v, dict):
                # Handle nested dictionaries
                json_results[k] = {
                    kk: float(vv) if isinstance(vv, (np.float32, np.float64)) else vv
                    for kk, vv in v.items()
                }
            else:
                json_results[k] = v
        
        json.dump(json_results, f, indent=2, default=str)
    
    if verbose:
        print(f"Saved results: {out/'results.json'}")
        print("\n=== Analysis Complete ===")


@main.command()
@click.option("--out", "outdir", type=click.Path(file_okay=False), default="screen_out",
              help="Output directory for data products and figures.")
@click.option("--skymap-dir", type=click.Path(file_okay=False), default=None,
              help="Sky-map download cache (default: <out>/skymaps).")
@click.option("--time-window", default="32 months ago .. now",
              help="GraceDB 'created:' window for the event query.")
@click.option("--far-threshold", type=float, default=2.3e-5,
              help="Keep events with false alarm rate below this (events/sec).")
@click.option("--bbh-threshold", type=float, default=0.1,
              help="Keep events with p_astro(BBH) above this. Use -1 to disable BBH filtering.")
@click.option("--workers", type=int, default=1,
              help="Worker processes for the all-pairs overlap computation.")
@click.option("--top-n", type=int, default=30, help="Number of ranked pairs to write to CSV.")
@click.option("--min-days-apart", type=int, default=0,
              help="Drop pairs fewer than N days apart (use 1 to exclude same-day duplicates).")
@click.option("--plot-pairs", type=int, default=1,
              help="Render joint sky-map plots for this many top-ranked pairs.")
@click.option("--reuse-overlaps", type=click.Path(exists=True, dir_okay=False), default=None,
              help="Skip query/download/compute and re-rank/re-plot from this overlaps JSON.")
def screen(outdir, skymap_dir, time_window, far_threshold, bbh_threshold,
           workers, top_n, min_days_apart, plot_pairs, reuse_overlaps):
    """Screen GW event pairs for high sky-map overlap (lensing candidates)."""
    # Imported here so the lightweight 'odds' path does not require the
    # screening extra (ligo.skymap, healpy, ligo-gracedb, ...).
    from gwassociation.screening.run import run_screening

    summary = run_screening(
        out_dir=outdir,
        skymap_dir=skymap_dir,
        time_window=time_window,
        far_threshold=far_threshold,
        bbh_threshold=None if bbh_threshold is not None and bbh_threshold < 0 else bbh_threshold,
        workers=workers,
        top_n=top_n,
        min_days_apart=min_days_apart,
        n_plot_pairs=plot_pairs,
        reuse_overlaps=reuse_overlaps,
    )

    top = summary["top_pairs"]
    click.echo("\nTop pairs by overlap:")
    for row in top[:10]:
        pv = row.get("pvalue")
        pv_str = f"  p={pv:.2e}" if pv is not None and pv == pv else ""
        click.echo(
            f"  {row['rank']:>2}. {row['event1']:<12} {row['event2']:<12} "
            f"O={row['overlap']:>10.3f}{pv_str}"
        )

    # List the files actually written into the output directory so the user
    # knows exactly what to open (reused inputs living elsewhere are excluded).
    out_path = pathlib.Path(outdir).resolve()
    outputs = summary.get("outputs", {})
    candidates = [outputs.get("top_pairs_csv"), outputs.get("overlaps_json")]
    candidates += list((outputs.get("figures") or {}).values())
    candidates.append(str(out_path / "summary.json"))
    written = []
    for candidate in candidates:
        if not candidate:
            continue
        path = pathlib.Path(candidate).resolve()
        if path.parent == out_path and path.exists() and path.name not in written:
            written.append(path.name)

    click.echo(f"\nData products written to {outdir}/:")
    for name in written:
        click.echo(f"  {name}")


@main.command(name="pair-odds")
@click.argument("skymap1", type=click.Path(exists=True, dir_okay=False))
@click.argument("skymap2", type=click.Path(exists=True, dir_okay=False))
@click.option("--model", type=click.Choice(["kilonova", "grb", "afterglow"]),
              default="kilonova", help="EM temporal model (used when times are given).")
@click.option("--gw-time", type=float, default=None,
              help="Event-1 GPS time; with --event2-time activates the temporal term.")
@click.option("--event2-time", type=float, default=None, help="Event-2 GPS time.")
@click.option("--no-distance", is_flag=True,
              help="Score sky position only (skip the distance term).")
@click.option("--out", "outfile", type=click.Path(dir_okay=False), default=None,
              help="Optional path to write the full result as JSON.")
def pair_odds(skymap1, skymap2, model, gw_time, event2_time, no_distance, outfile):
    """Point-source association odds for a PAIR of GW events.

    Reduces SKYMAP2 to a point source (its peak position + distance-derived
    redshift) and scores it against SKYMAP1 with the same odds path used for
    GW170817. Prints the I_omega / I_dl / I_t / odds / P breakdown.
    """
    from gwassociation.screening.bridge import pair_point_odds

    name1 = pathlib.Path(skymap1).name.split(".")[0]
    name2 = pathlib.Path(skymap2).name.split(".")[0]

    r = pair_point_odds(
        skymap1, skymap2, name1=name1, name2=name2, em_model=model,
        gw_time=gw_time, event2_time=event2_time, include_distance=not no_distance,
    )
    t = r["_transient"]

    click.echo(f"\n=== {name1}  vs  {name2} (as point source) ===")
    z_str = f", z={t['z']:.4f}" if "z" in t else " (distance off)"
    click.echo(f"  point: RA={t['ra']:.2f}, Dec={t['dec']:.2f}{z_str}")
    click.echo(f"  I_omega        : {r['I_omega']:.3e}")
    click.echo(f"  I_dl           : {r['I_dl']:.3e}")
    click.echo(f"  I_t            : {r['I_t']:.3e}")
    click.echo(f"  Bayes factor   : {r['bayes_factor']:.3e}")
    click.echo(f"  posterior odds : {r['posterior_odds']:.3e}")
    click.echo(f"  log10 odds     : {r['log_posterior_odds']:.2f}")
    click.echo(f"  P(assoc)       : {r['confidence']:.1%}")

    if outfile:
        payload = {k: v for k, v in r.items() if not k.startswith("_")}
        payload["transient"] = t
        with open(outfile, "w") as handle:
            json.dump(payload, handle, indent=2, default=str)
        click.echo(f"\nWrote {outfile}")


if __name__ == "__main__":
    main()