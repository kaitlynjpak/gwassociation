import os
import json
import pathlib
import click
import numpy as np

@click.command()
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
def main(gw_file, secondary_skymap, ra, dec, z, z_err, ttime,
         secondary_time, gw_time, model, outdir, verbose):
    """Run GW-EM association analysis"""
    
    # Only import here to avoid issues if package not fully installed
    from gw_assoc import Association
    from gw_assoc.plots import plot_association_summary
    
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

if __name__ == "__main__":
    main()