#!/usr/bin/env python3
"""
Test the GW-EM association framework with real GW data
"""

import sys
from pathlib import Path

# Make sure the package is in the path
sys.path.insert(0, 'src')

from gw_assoc import Association

def test_with_real_skymap():
    """
    Test with a real GW skymap file
    
    Update the skymap filename and transient parameters 
    with your actual data
    """
    
    # ============================================================
    # OPTION 1: If you downloaded a real skymap
    # ============================================================
    
    # Uncomment and modify this section if you have a real skymap:
    skymap_file = "S190521_bayestar.fits.gz,0"  # Replace with your file
    
    if not Path(skymap_file).exists():
        print(f"✗ Skymap file not found: {skymap_file}")
        print("\nTo get a real skymap:")
        print("1. Go to: https://gracedb.ligo.org/superevents/public/O3/")
        print("2. Click on an event (e.g., S190425z)")
        print("3. Download the bayestar.fits.gz file")
        print("4. Put it in this directory and update the filename above")
        return
    
    print(f"✓ Using skymap: {skymap_file}")
    
    # Create association with real data
    assoc = Association(
        gw_file=skymap_file,
        transient_info={
            "name": "AT2019ebq",  # Your transient name
            "ra": 192.42625,      # Replace with your transient's RA
            "dec": 34.82472,     # Replace with your transient's Dec
            "z": 0.438,       # Replace with your transient's redshift
            "z_err": None,   # Redshift uncertainty
            "time": None,    # Your transient detection time
            "gw_time": 1242442967.447, # GW trigger time from GraceDB
            "magnitude": None,
        }
    )
    
    # ============================================================
    # OPTION 2: Use the test skymap (default)
    # ============================================================
    
    skymap_file = "test_data/test_skymap.fits"
    
    if not Path(skymap_file).exists():
        print("✗ Test skymap not found!")
        print("  Run 'python test_gw_assoc.py' first to create test data")
        return
    
    print(f"✓ Using test skymap: {skymap_file}")
    print("  (To use real data, see instructions in the code)\n")
    
    # Create association with test data
    # This transient is well-matched to the test skymap
    assoc = Association(
        gw_file=skymap_file,
        transient_info={
            "ra": 120.5,      # Near center of test skymap
            "dec": -30.0,     # Near center of test skymap
            "z": 0.050,       # Reasonable redshift
            "z_err": 0.005,   
            "time": 1234567890.0,    
            "gw_time": 1234567880.0,  # 10 seconds before transient
            "magnitude": 18.5,
            "name": "AT2024test"
        }
    )
    
    # ============================================================
    # Run the analysis
    # ============================================================
    
    print("=" * 60)
    print("GW-EM ASSOCIATION ANALYSIS")
    print("=" * 60)
    
    # Compute association with kilonova model
    results = assoc.compute_odds(
        em_model='kilonova',
        chance_coincidence_rate=1e-5,
        H0_uncertainty=7.0
    )
    
    # Display results
    print(f"\nTransient: {assoc.transient.name or 'Unknown'}")
    print(f"Position:  RA = {assoc.transient.ra:.3f}°, Dec = {assoc.transient.dec:.3f}°")
    if assoc.transient.z:
        print(f"Redshift:  z = {assoc.transient.z:.4f} ± {assoc.transient.z_err or 0:.4f}")
        from astropy.cosmology import Planck15
        import astropy.units as u
        d_L = Planck15.luminosity_distance(assoc.transient.z).to(u.Mpc).value
        print(f"Distance:  {d_L:.1f} Mpc")
    
    print(f"\nGW Event:")
    print(f"Time:      {assoc.gw.event_time}")
    print(f"Time delay: {(assoc.transient.time - assoc.gw.event_time):.1f} seconds")
    
    print(f"\n" + "-" * 40)
    print("OVERLAP INTEGRALS:")
    print("-" * 40)
    print(f"Spatial (I_Ω):   {results['I_omega']:.3e}  {'✓' if results['I_omega'] > 1e-4 else '✗'}")
    print(f"Distance (I_DL):  {results['I_dl']:.3e}  {'✓' if results['I_dl'] > 0.1 else '✗'}")
    print(f"Temporal (I_t):   {results['I_t']:.3e}  {'✓' if results['I_t'] > 0.01 else '✗'}")
    
    print(f"\n" + "-" * 40)
    print("STATISTICAL RESULTS:")
    print("-" * 40)
    print(f"Bayes Factor:     {results['bayes_factor']:.3e}")
    print(f"Posterior Odds:   {results['posterior_odds']:.3e}")
    print(f"Log₁₀ Odds:       {results['log_posterior_odds']:.2f}")
    
    print(f"\n" + "-" * 40)
    print("ASSOCIATION PROBABILITY:")
    print("-" * 40)
    prob = results['confidence']
    print(f"P(Associated) = {prob:.1%}")
    
    # Interpretation
    if prob > 0.95:
        print("\n🌟 STRONG ASSOCIATION - Very likely the EM counterpart!")
    elif prob > 0.5:
        print("\n✓ PROBABLE ASSOCIATION - Good candidate for follow-up")
    elif prob > 0.1:
        print("\n⚠ WEAK ASSOCIATION - Possible but unlikely")
    else:
        print("\n✗ NO ASSOCIATION - Very unlikely to be related")
    
    # Generate plots
    print(f"\n" + "-" * 40)
    print("GENERATING PLOTS:")
    print("-" * 40)
    
    try:
        output_dir = Path("gw_analysis_output")
        output_dir.mkdir(exist_ok=True)
        
        # Skymap plot
        skymap_file = output_dir / "skymap.png"
        assoc.plot_skymap(str(skymap_file))
        print(f"✓ Skymap saved to: {skymap_file}")
        
        # Summary plot
        from gw_assoc.plots import plot_association_summary
        summary_file = output_dir / "summary.png"
        plot_association_summary(results, str(summary_file))
        print(f"✓ Summary saved to: {summary_file}")
        
    except Exception as e:
        print(f"⚠ Could not generate plots: {e}")
    
    # Save results
    import json
    results_file = output_dir / "results.json"
    with open(results_file, 'w') as f:
        # Convert numpy types for JSON
        json_results = {k: float(v) if isinstance(v, (float, int)) and v != v else v 
                       for k, v in results.items()}
        json.dump(json_results, f, indent=2, default=str)
    print(f"✓ Results saved to: {results_file}")
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)

def test_multiple_transients():
    """
    Test with multiple transient candidates
    """
    print("\n" + "=" * 60)
    print("TESTING MULTIPLE CANDIDATES")
    print("=" * 60)
    
    # Define your candidates
    candidates = [
        {
            "name": "ZTF19aadyppr",
            "ra": 120.1,
            "dec": -30.2,
            "z": 0.048,
            "z_err": 0.003,
            "time": 1234567885.0,
            "magnitude": 17.8
        },
        {
            "name": "ATLAS19dfq",
            "ra": 121.5,
            "dec": -29.5,
            "z": 0.052,
            "z_err": 0.005,
            "time": 1234567890.0,
            "magnitude": 18.2
        },
        {
            "name": "MASTER_OT_J1234",
            "ra": 119.8,
            "dec": -31.0,
            "z": 0.055,
            "z_err": 0.004,
            "time": 1234567900.0,
            "magnitude": 18.5
        }
    ]
    
    # Use test skymap
    skymap_file = "test_data/test_skymap.fits"
    if not Path(skymap_file).exists():
        print("✗ Test skymap not found! Run test_gw_assoc.py first")
        return
    
    # Create association
    assoc = Association(
        gw_file=skymap_file,
        transient_info={"ra": 0.0, "dec": 0.0, "gw_time": 1234567880.0}
    )

    
    # Rank candidates
    print("\nRanking candidates...")
    rankings = assoc.rank_candidates(candidates)
    
    # Display results
    print("\n" + "-" * 70)
    print(f"{'Rank':<6} {'Name':<20} {'RA':<8} {'Dec':<8} {'z':<8} {'P(Assoc)':<12}")
    print("-" * 70)
    
    for i, r in enumerate(rankings, 1):
        c = r['candidate']
        print(f"{i:<6} {c.name:<20} {c.ra:<8.2f} {c.dec:<8.2f} "
              f"{c.z:<8.4f} {r['probability']:<12.1%}")
    
    print("\n✓ Best candidate: " + rankings[0]['candidate'].name)

if __name__ == "__main__":
    # Run the single transient test
    test_with_real_skymap()
    
    # Optionally test multiple candidates
    print("\n" + "=" * 60)
    choice = input("Test multiple candidates? (y/n): ")
    if choice.lower() == 'y':
        test_multiple_transients()