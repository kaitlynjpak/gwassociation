"""
Example usage of the GW-EM Association Framework
This demonstrates how to use the framework with real-world scenarios
"""

import numpy as np
from pathlib import Path

# ============================================================================
# EXAMPLE 1: Simple Single Transient Analysis
# ============================================================================

def example_single_transient():
    """Analyze a single transient candidate"""
    print("\n" + "="*60)
    print("EXAMPLE 1: Single Transient Analysis")
    print("="*60)
    
    from gw_assoc import Association
    
    # Create association with GW event and transient
    assoc = Association(
        gw_file="path/to/S190425z_skymap.fits",  # Replace with actual skymap
        transient_info={
            "ra": 120.5,      # Right ascension in degrees
            "dec": -30.2,     # Declination in degrees  
            "z": 0.048,       # Redshift
            "z_err": 0.003,   # Redshift uncertainty
            "time": 1240215505.0,  # Detection time (GPS)
            "gw_time": 1240215503.0,  # GW event time (GPS)
            "magnitude": 17.8,
            "filter_band": "r",
            "name": "AT2019test"
        }
    )
    
    # Compute association odds
    results = assoc.compute_odds(
        em_model='kilonova',  # or 'grb', 'afterglow'
        chance_coincidence_rate=1e-5,
        H0_uncertainty=7.0  # km/s/Mpc
    )
    
    # Print results
    print(f"\nTransient: AT2019test")
    print(f"Position: RA={120.5}°, Dec={-30.2}°")
    print(f"Redshift: z={0.048} ± 0.003")
    print(f"\nOverlap Integrals:")
    print(f"  Spatial (I_Ω):  {results['I_omega']:.3e}")
    print(f"  Distance (I_DL): {results['I_dl']:.3e}")
    print(f"  Temporal (I_t):  {results['I_t']:.3e}")
    print(f"\nStatistics:")
    print(f"  Bayes Factor:    {results['bayes_factor']:.3e}")
    print(f"  Posterior Odds:  {results['posterior_odds']:.3e}")
    print(f"  P(Associated):   {results['confidence']:.1%}")
    print(f"\nDecision: {'✓ ASSOCIATED' if results['associated'] else '✗ NOT ASSOCIATED'}")
    
    # Generate plots
    assoc.plot_skymap("example_skymap.png")
    print(f"\nSaved skymap to: example_skymap.png")

# ============================================================================
# EXAMPLE 2: Multiple Candidates Ranking
# ============================================================================

def example_multiple_candidates():
    """Rank multiple transient candidates"""
    print("\n" + "="*60)
    print("EXAMPLE 2: Ranking Multiple Candidates")
    print("="*60)
    
    from gw_assoc import Association
    
    # Define multiple candidates
    candidates = [
        {
            "name": "AT2019aaa",
            "ra": 120.5,
            "dec": -30.2,
            "z": 0.048,
            "z_err": 0.003,
            "time": 1240215505.0,
            "magnitude": 17.8
        },
        {
            "name": "AT2019bbb", 
            "ra": 121.2,
            "dec": -29.8,
            "z": 0.051,
            "z_err": 0.005,
            "time": 1240215600.0,
            "magnitude": 18.2
        },
        {
            "name": "AT2019ccc",
            "ra": 119.8,
            "dec": -31.0,
            "z": 0.12,  # Wrong distance
            "z_err": 0.01,
            "time": 1240215510.0,
            "magnitude": 19.5
        },
        {
            "name": "AT2019ddd",
            "ra": 200.0,  # Far from localization
            "dec": 45.0,
            "z": 0.05,
            "z_err": 0.003,
            "time": 1240215504.0,
            "magnitude": 18.0
        }
    ]
    
    # Create association
    assoc = Association(
        gw_file="path/to/skymap.fits",
        transient_info={"gw_time": 1240215503.0}
    )
    
    # Rank candidates
    rankings = assoc.rank_candidates(candidates)
    
    # Display rankings
    print("\nCandidate Rankings:")
    print("-" * 60)
    print(f"{'Rank':<6} {'Name':<12} {'RA':<8} {'Dec':<8} {'z':<8} {'P(Assoc)':<10} {'Decision'}")
    print("-" * 60)
    
    for i, r in enumerate(rankings):
        candidate = r['candidate']
        prob = r['probability']
        decision = "✓" if prob > 0.5 else "✗"
        
        print(f"{i+1:<6} {candidate.name:<12} "
              f"{candidate.ra:<8.2f} {candidate.dec:<8.2f} "
              f"{candidate.z:<8.3f} {prob:<10.1%} {decision}")
    
    # Find best candidate
    best = rankings[0]
    print(f"\nBest candidate: {best['candidate'].name} "
          f"with P(Associated) = {best['probability']:.1%}")

# ============================================================================
# EXAMPLE 3: Different EM Models
# ============================================================================

def example_different_models():
    """Compare different EM counterpart models"""
    print("\n" + "="*60)
    print("EXAMPLE 3: Comparing EM Models")
    print("="*60)
    
    from gw_assoc import Association
    
    # Same transient, different models
    transient_info = {
        "ra": 120.5,
        "dec": -30.2,
        "z": 0.05,
        "time": 1240215504.7,  # 1.7 seconds after GW
        "gw_time": 1240215503.0
    }
    
    assoc = Association("path/to/skymap.fits", transient_info)
    
    models = ['grb', 'kilonova', 'afterglow']
    
    print(f"\nTransient detected {1.7} seconds after GW")
    print("Model comparison:")
    print("-" * 40)
    
    for model in models:
        results = assoc.compute_odds(em_model=model)
        print(f"{model:10} → P(Associated) = {results['confidence']:.1%}")
    
    print("\nInterpretation:")
    print("- GRB: Best for prompt emission (seconds)")
    print("- Kilonova: Best for optical/NIR (hours-days)")
    print("- Afterglow: Best for late-time (days-weeks)")

# ============================================================================
# EXAMPLE 4: Real GW Events
# ============================================================================

def example_real_events():
    """Examples based on real GW events"""
    print("\n" + "="*60)
    print("EXAMPLE 4: Real GW Event Scenarios")
    print("="*60)
    
    # GW170817-like event (golden standard)
    print("\n1. GW170817-like (BNS at 40 Mpc):")
    print("   - Small sky area: ~30 deg²")
    print("   - Well-constrained distance: 40±8 Mpc")
    print("   - Kilonova detected after 11 hours")
    print("   Expected: I_Ω~0.01, I_DL~0.8, I_t~0.9")
    print("   → Very high association probability")
    
    # O4-like event (typical)
    print("\n2. Typical O4 BNS (200 Mpc):")
    print("   - Large sky area: ~500 deg²")
    print("   - Moderate distance constraint: 200±50 Mpc")
    print("   - Multiple candidates expected")
    print("   Expected: I_Ω~10⁻⁴, I_DL~0.3, I_t~0.5")
    print("   → Need careful ranking of candidates")
    
    # NSBH event
    print("\n3. NS-BH merger:")
    print("   - Very large sky area: ~1000 deg²")
    print("   - Poor distance constraint")
    print("   - Fainter EM counterpart expected")
    print("   → Challenging association")

# ============================================================================
# EXAMPLE 5: Workflow for Follow-up Campaign  
# ============================================================================

def example_followup_workflow():
    """Complete workflow for GW follow-up"""
    print("\n" + "="*60)
    print("EXAMPLE 5: GW Follow-up Workflow")
    print("="*60)
    
    print("""
    STEP 1: GW ALERT RECEIVED
    ├─ Download skymap from GraceDB
    ├─ Extract event parameters (time, distance, localization)
    └─ Estimate search area (90% credible region)
    
    STEP 2: TELESCOPE OBSERVATIONS
    ├─ Calculate tiles needed for coverage
    ├─ Prioritize by probability and galaxy density
    └─ Begin observations (typically within hours)
    
    STEP 3: TRANSIENT DISCOVERY
    ├─ Image subtraction to find new sources
    ├─ Photometry and preliminary classification
    └─ Compile candidate list
    
    STEP 4: ASSOCIATION ANALYSIS (this framework)
    ├─ Load GW skymap
    ├─ For each candidate:
    │   ├─ Calculate I_Ω (spatial overlap)
    │   ├─ Calculate I_DL (distance overlap) 
    │   ├─ Calculate I_t (temporal overlap)
    │   └─ Compute posterior odds
    ├─ Rank all candidates
    └─ Identify most likely counterpart
    
    STEP 5: FOLLOW-UP DECISIONS
    ├─ Spectroscopy for top candidates
    ├─ Multi-wavelength observations
    └─ Report to GCN/TNS
    """)
    
    print("\nExample code for Step 4:")
    print("-" * 40)
    print("""
    # Load skymap and candidates
    from gw_assoc import Association
    import json
    
    # Read candidate list from observers
    with open('candidates.json') as f:
        candidates = json.load(f)
    
    # Create association and rank
    assoc = Association('S230518h_skymap.fits', 
                       {'gw_time': gw_trigger_time})
    rankings = assoc.rank_candidates(candidates)
    
    # Get top 3 for spectroscopy
    top_targets = rankings[:3]
    for target in top_targets:
        print(f"{target['candidate'].name}: "
              f"P={target['probability']:.1%}")
    """)

# ============================================================================
# EXAMPLE 6: Advanced Features
# ============================================================================

def example_advanced_features():
    """Demonstrate advanced features of the framework"""
    print("\n" + "="*60)
    print("EXAMPLE 6: Advanced Features")
    print("="*60)
    
    from gw_assoc import Association
    from gw_assoc.analysis.los import DistanceOverlap
    from gw_assoc.density import false_alarm_probability
    from gw_assoc.fov import tiles_needed, observation_time
    
    print("\n1. DISTANCE WITH PECULIAR VELOCITIES")
    print("-" * 40)
    
    # Account for peculiar velocity uncertainty
    distance_calc = DistanceOverlap()
    print("For nearby galaxy (z=0.01):")
    print("  Peculiar velocity: ~300 km/s")
    print("  Fractional error: ~10%")
    print("  This significantly affects I_DL!")
    
    print("\n2. FALSE ALARM CALCULATION")
    print("-" * 40)
    
    # Calculate false alarm probability
    n_candidates = 5
    search_area = 500  # deg²
    search_time = 3  # days
    
    p_false = false_alarm_probability(n_candidates, search_area, search_time)
    print(f"Found {n_candidates} candidates in {search_area} deg²")
    print(f"False alarm probability: {p_false:.2%}")
    
    print("\n3. OBSERVATION PLANNING")
    print("-" * 40)
    
    skymap_area = 500  # deg²
    
    for telescope in ['ZTF', 'DECam', 'LSST']:
        n_tiles = tiles_needed(skymap_area, telescope)
        obs_time = observation_time(skymap_area, telescope)
        print(f"{telescope:10} → {n_tiles:3d} tiles, {obs_time:.1f} hours")
    
    print("\n4. CUSTOM PRIORS")
    print("-" * 40)
    
    # Use different priors based on GW source type
    print("Prior odds by source type:")
    print("  BNS (NS-NS):     1.0  (EM expected)")
    print("  NSBH (NS-BH):    0.1  (EM possible)")  
    print("  BBH (BH-BH):     0.01 (EM unlikely)")
    
    # Example with NSBH prior
    assoc = Association("skymap.fits", {"ra": 120, "dec": -30, "z": 0.05})
    results_nsbh = assoc.compute_odds(prior_odds=0.1)
    print(f"\nWith NSBH prior: P(Associated) = {results_nsbh['confidence']:.1%}")

# ============================================================================
# MAIN: Run Examples
# ============================================================================

def main():
    """Run all examples"""
    import sys
    
    print("""
╔══════════════════════════════════════════════════════════╗
║          GW-EM ASSOCIATION FRAMEWORK EXAMPLES           ║
╚══════════════════════════════════════════════════════════╝
    """)
    
    print("This script demonstrates how to use the framework.")
    print("Note: Some examples need actual skymap files to run.\n")
    
    # Check if we can import the package
    try:
        import gw_assoc
        print("✓ GW-Assoc package found\n")
    except ImportError:
        print("✗ Cannot import gw_assoc package!")
        print("  Make sure it's in your Python path:")
        print("  export PYTHONPATH=$PYTHONPATH:/path/to/gw_assoc")
        sys.exit(1)
    
    # Menu
    examples = {
        '1': ('Single Transient Analysis', example_single_transient),
        '2': ('Multiple Candidates Ranking', example_multiple_candidates),
        '3': ('Different EM Models', example_different_models),
        '4': ('Real GW Event Scenarios', example_real_events),
        '5': ('Follow-up Workflow', example_followup_workflow),
        '6': ('Advanced Features', example_advanced_features),
    }
    
    print("Available examples:")
    for key, (name, _) in examples.items():
        print(f"  {key}. {name}")
    print("  A. Run all examples")
    print("  Q. Quit")
    
    while True:
        choice = input("\nSelect an example (1-6, A, Q): ").strip().upper()
        
        if choice == 'Q':
            break
        elif choice == 'A':
            for name, func in examples.values():
                try:
                    func()
                except Exception as e:
                    print(f"\n⚠ Example '{name}' failed: {e}")
                    print("  (This may be due to missing skymap files)")
            break
        elif choice in examples:
            try:
                examples[choice][1]()
            except Exception as e:
                print(f"\n⚠ Example failed: {e}")
                print("  (This may be due to missing skymap files)")
        else:
            print("Invalid choice. Please select 1-6, A, or Q.")
    
    print("\n" + "="*60)
    print("Examples completed!")
    print("\nNext steps:")
    print("1. Download real GW skymaps from: https://gracedb.ligo.org/superevents/")
    print("2. Prepare your transient candidate list")
    print("3. Run the association analysis")
    print("4. Rank candidates and plan follow-up observations")

if __name__ == "__main__":
    main()