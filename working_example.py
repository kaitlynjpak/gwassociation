import os
import sys
from pathlib import Path

def run_working_example():
    """Run a complete working example with test data"""
    
    print("=" * 60)
    print("GW-EM ASSOCIATION WORKING EXAMPLE")
    print("=" * 60)
    
    # Check if test data exists
    if not Path("test_data/test_skymap.fits").exists():
        print("\n✗ Test data not found!")
        print("Please run: python test_gw_assoc.py")
        print("This will create the test data needed for this example.")
        return False
    
    print("\n✓ Test data found. Running example...\n")
    
    # Import the framework
    try:
        from gw_assoc import Association
        from gw_assoc.plots import plot_association_summary
    except ImportError as e:
        print(f"✗ Import error: {e}")
        print("Make sure gw_assoc is in your Python path")
        return False
    
    # ============================================================
    # EXAMPLE 1: Single Transient (Good Match)
    # ============================================================
    print("EXAMPLE 1: Good Match (close to GW localization)")
    print("-" * 40)
    
    assoc1 = Association(
        gw_file="test_data/test_skymap.fits",
        transient_info={
            "ra": 120.5,      # Close to skymap center (120°)
            "dec": -29.8,     # Close to skymap center (-30°)
            "z": 0.05,
            "z_err": 0.005,
            "time": 1234567891.0,    # 11 seconds after GW
            "gw_time": 1234567880.0,
            "magnitude": 17.8,
            "filter_band": "r",
            "name": "AT2024good"
        }
    )
    
    # Compute odds with kilonova model
    results1 = assoc1.compute_odds(
        em_model='kilonova',
        chance_coincidence_rate=1e-5
    )
    
    print(f"Transient: AT2024good")
    print(f"Position: RA=120.5°, Dec=-29.8°")
    print(f"Redshift: z=0.050 ± 0.005")
    print(f"Time delay: 11 seconds")
    print(f"\nResults:")
    print(f"  Spatial overlap (I_Ω):  {results1['I_omega']:.3e}")
    print(f"  Distance overlap (I_DL): {results1['I_dl']:.3e}")
    print(f"  Temporal overlap (I_t):  {results1['I_t']:.3e}")
    print(f"  Bayes Factor:           {results1['bayes_factor']:.3e}")
    print(f"  P(Associated):          {results1['confidence']:.1%}")
    print(f"  Decision:               {'✓ ASSOCIATED' if results1['associated'] else '✗ NOT ASSOCIATED'}")
    
    # Generate skymap plot
    try:
        assoc1.plot_skymap("example1_skymap.png")
        print(f"\n✓ Saved skymap plot: example1_skymap.png")
    except Exception as e:
        print(f"\n⚠ Could not create skymap plot: {e}")
    
    # ============================================================
    # EXAMPLE 2: Poor Match (far from localization)
    # ============================================================
    print("\n" + "=" * 60)
    print("EXAMPLE 2: Poor Match (far from GW localization)")
    print("-" * 40)
    
    assoc2 = Association(
        gw_file="test_data/test_skymap.fits",
        transient_info={
            "ra": 200.0,      # Far from skymap center
            "dec": 45.0,      # Far from skymap center
            "z": 0.08,        # Wrong distance
            "z_err": 0.01,
            "time": 1234568000.0,    # 120 seconds after GW
            "gw_time": 1234567880.0,
            "magnitude": 19.2,
            "filter_band": "g",
            "name": "AT2024far"
        }
    )
    
    results2 = assoc2.compute_odds(em_model='kilonova')
    
    print(f"Transient: AT2024far")
    print(f"Position: RA=200.0°, Dec=45.0°")
    print(f"Redshift: z=0.080 ± 0.010")
    print(f"Time delay: 120 seconds")
    print(f"\nResults:")
    print(f"  Spatial overlap (I_Ω):  {results2['I_omega']:.3e}")
    print(f"  Distance overlap (I_DL): {results2['I_dl']:.3e}")
    print(f"  Temporal overlap (I_t):  {results2['I_t']:.3e}")
    print(f"  P(Associated):          {results2['confidence']:.1%}")
    print(f"  Decision:               {'✓ ASSOCIATED' if results2['associated'] else '✗ NOT ASSOCIATED'}")
    
    # ============================================================
    # EXAMPLE 3: Multiple Candidates Ranking
    # ============================================================
    print("\n" + "=" * 60)
    print("EXAMPLE 3: Ranking Multiple Candidates")
    print("-" * 40)
    
    candidates = [
        {
            "name": "AT2024best",
            "ra": 120.2,
            "dec": -30.1,
            "z": 0.049,
            "z_err": 0.003,
            "time": 1234567890.0,
            "magnitude": 17.5
        },
        {
            "name": "AT2024okay",
            "ra": 121.5,
            "dec": -29.0,
            "z": 0.052,
            "z_err": 0.005,
            "time": 1234567900.0,
            "magnitude": 18.0
        },
        {
            "name": "AT2024poor",
            "ra": 125.0,
            "dec": -35.0,
            "z": 0.045,
            "z_err": 0.004,
            "time": 1234568000.0,
            "magnitude": 18.5
        },
        {
            "name": "AT2024bad",
            "ra": 180.0,
            "dec": 0.0,
            "z": 0.1,
            "z_err": 0.01,
            "time": 1234569000.0,
            "magnitude": 19.0
        }
    ]
    
    # Add gw_time for ranking
    for c in candidates:
        c['gw_time'] = 1234567880.0
    
    # Create association and rank
    assoc3 = Association(
        gw_file="test_data/test_skymap.fits",
        transient_info={"gw_time": 1234567880.0}
    )
    
    rankings = assoc3.rank_candidates(candidates)
    
    print("\nCandidate Rankings:")
    print(f"{'Rank':<6} {'Name':<12} {'RA':<8} {'Dec':<8} {'z':<8} {'P(Assoc)':<12} {'Decision'}")
    print("-" * 70)
    
    for i, r in enumerate(rankings):
        candidate = r['candidate']
        prob = r['probability']
        decision = "✓ ASSOCIATED" if prob > 0.5 else "✗ UNLIKELY"
        
        print(f"{i+1:<6} {candidate.name:<12} "
              f"{candidate.ra:<8.1f} {candidate.dec:<8.1f} "
              f"{candidate.z:<8.3f} {prob:<12.1%} {decision}")
    
    best = rankings[0]
    print(f"\nBest candidate: {best['candidate'].name} with P = {best['probability']:.1%}")
    
    # Try to create ranking plot
    try:
        from gw_assoc.plots import plot_candidate_ranking
        plot_candidate_ranking(rankings, "example_ranking.png")
        print("✓ Saved ranking plot: example_ranking.png")
    except Exception as e:
        print(f"⚠ Could not create ranking plot: {e}")
    
    # ============================================================
    # EXAMPLE 4: Different EM Models
    # ============================================================
    print("\n" + "=" * 60)
    print("EXAMPLE 4: Comparing Different EM Models")
    print("-" * 40)
    
    # Transient detected quickly after GW
    fast_transient = {
        "ra": 120.5,
        "dec": -30.0,
        "z": 0.05,
        "time": 1234567882.0,  # 2 seconds after GW
        "gw_time": 1234567880.0
    }
    
    print("Transient detected 2 seconds after GW")
    print("Comparing models:")
    print("-" * 30)
    
    models = ['grb', 'kilonova', 'afterglow']
    for model in models:
        assoc_model = Association("test_data/test_skymap.fits", fast_transient.copy())
        results = assoc_model.compute_odds(em_model=model)
        print(f"{model:10} → I_t = {results['I_t']:.3f}, P(Assoc) = {results['confidence']:.1%}")
    
    print("\nInterpretation:")
    print("- GRB model favors prompt detection (seconds)")
    print("- Kilonova model expects delay (hours-days)")
    print("- This rapid detection suggests GRB/prompt emission")
    
    # ============================================================
    # SUMMARY
    # ============================================================
    print("\n" + "=" * 60)
    print("EXAMPLE COMPLETE!")
    print("=" * 60)
    print("\nKey Takeaways:")
    print("1. Position match (I_Ω) is crucial - being in the skymap matters most")
    print("2. Distance agreement (I_DL) helps confirm association")
    print("3. Time delay (I_t) depends on the EM model")
    print("4. Multiple candidates should be ranked by posterior odds")
    print("\nGenerated files:")
    print("- example1_skymap.png (if healpy available)")
    print("- example_ranking.png (if plotting works)")
    
    return True

if __name__ == "__main__":
    # Run the working example
    success = run_working_example()
    
    if not success:
        sys.exit(1)
    
    print("\n✓ All examples completed successfully!")
    print("\nNext steps:")
    print("1. Try with real GW skymaps from GraceDB")
    print("2. Use your own transient candidates")
    print("3. Adjust the EM models and priors for your science case")