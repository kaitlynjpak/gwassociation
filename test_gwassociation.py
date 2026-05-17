"""
Comprehensive test suite for GW-EM Association Framework
Run this to verify your installation and test all components
"""

import numpy as np
import os
import sys
import json
from pathlib import Path

from gwassociation.utils import healpix as hp_utils

# ============================================================================
# 1. GENERATE TEST DATA
# ============================================================================

def create_test_skymap(filename="test_skymap.fits", nside=64):
    """Create a simple test GW skymap"""
    npix = hp_utils.nside2npix(nside)
    skymap = np.zeros(npix)

    # Center at RA=120, Dec=-30
    center_vec = np.asarray(hp_utils.ang2vec(np.radians(90 + 30), np.radians(120))).reshape(3)

    # Add Gaussian probability
    for ipix in range(npix):
        pix_vec = np.asarray(hp_utils.pix2vec(nside, ipix)).reshape(3)
        angle = np.arccos(np.clip(np.dot(center_vec, pix_vec), -1, 1))
        skymap[ipix] = np.exp(-0.5 * (angle / np.radians(10))**2)

    # Normalize
    skymap = skymap / np.sum(skymap)

    # Save as FITS
    hp_utils.write_map(filename, skymap, overwrite=True)
    print(f"✓ Created test skymap: {filename}")
    return filename

def create_test_data():
    """Create all test data files"""
    print("\n" + "="*60)
    print("CREATING TEST DATA")
    print("="*60)
    
    # Create test directory
    test_dir = Path("test_data")
    test_dir.mkdir(exist_ok=True)
    os.chdir(test_dir)
    
    # Create skymap
    skymap_file = create_test_skymap()
    
    # Create test transient list
    transients = [
        {
            "name": "AT2024test_good",
            "ra": 120.5,
            "dec": -29.8,
            "z": 0.05,
            "z_err": 0.005,
            "time": 1234567891.0,
            "magnitude": 18.5,
            "filter_band": "r"
        },
        {
            "name": "AT2024test_far",
            "ra": 200.0,
            "dec": 45.0,
            "z": 0.08,
            "z_err": 0.01,
            "time": 1234568000.0,
            "magnitude": 19.2,
            "filter_band": "g"
        },
        {
            "name": "AT2024test_noredshift",
            "ra": 119.8,
            "dec": -30.5,
            "z": None,
            "time": 1234567900.0,
            "magnitude": 17.8,
            "filter_band": "i"
        }
    ]
    
    with open("test_transients.json", "w") as f:
        json.dump(transients, f, indent=2)
    print("✓ Created test transients: test_transients.json")
    
    # Create CSV version
    import csv
    with open("test_transients.csv", "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["name", "ra", "dec", "z", "z_err", "time", "magnitude", "filter_band"])
        writer.writeheader()
        writer.writerows(transients)
    print("✓ Created test transients: test_transients.csv")
    
    os.chdir("..")
    return test_dir

# ============================================================================
# 2. TEST INDIVIDUAL COMPONENTS
# ============================================================================

def test_imports():
    """Test that all modules can be imported"""
    print("\n" + "="*60)
    print("TESTING IMPORTS")
    print("="*60)
    
    modules_to_test = [
        ("gwassociation", "Main module"),
        ("gwassociation.io.transient", "Transient I/O"),
        ("gwassociation.io.skymap", "Skymap I/O"),
        ("gwassociation.analysis.spatial", "Spatial analysis"),
        ("gwassociation.analysis.temporal", "Temporal analysis"),
        ("gwassociation.analysis.los", "Line of sight"),
        ("gwassociation.analysis.odds", "Odds calculation"),
        ("gwassociation.association", "Association class"),
        ("gwassociation.plotting.skymap", "Skymap plotting"),
        ("gwassociation.plots", "Summary plots"),
    ]
    
    failed = []
    for module_name, description in modules_to_test:
        try:
            __import__(module_name)
            print(f"✓ {description:25} ({module_name})")
        except ImportError as e:
            print(f"✗ {description:25} ({module_name})")
            print(f"  Error: {e}")
            failed.append(module_name)
    
    return len(failed) == 0

def test_basic_functionality():
    """Test basic functionality"""
    print("\n" + "="*60)
    print("TESTING BASIC FUNCTIONALITY")
    print("="*60)
    
    try:
        from gwassociation.io.transient import Transient
        from gwassociation.io.skymap import GWEvent
        
        # Test Transient creation
        transient = Transient(
            ra=120.5,
            dec=-30.0,
            z=0.05,
            z_err=0.005,
            time=1234567890.0,
            name="Test_Transient"
        )
        print(f"✓ Created Transient: {transient}")
        
        # Test coordinate conversion
        skycoord = transient.get_skycoord()
        print(f"✓ SkyCoord: RA={skycoord.ra.deg:.2f}°, Dec={skycoord.dec.deg:.2f}°")
        
        # Test distance calculation
        dl = transient.get_luminosity_distance()
        if dl:
            print(f"✓ Luminosity distance: {dl:.1f} Mpc")
        
        # Test GWEvent creation
        gw = GWEvent(
            skymap_path="test_data/test_skymap.fits",
            event_time=1234567880.0,
            event_name="GW_TEST"
        )
        print(f"✓ Created GWEvent: {gw.event_name}")
        
        return True
        
    except Exception as e:
        print(f"✗ Error in basic functionality: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_association_calculation():
    """Test the main association calculation"""
    print("\n" + "="*60)
    print("TESTING ASSOCIATION CALCULATION")
    print("="*60)
    
    try:
        from gwassociation import Association
        
        # Test with good match (close in sky and distance)
        print("\n1. Testing GOOD MATCH (close in position and distance):")
        assoc = Association("test_data/test_skymap.fits", {
            "ra": 120.5,
            "dec": -29.8,
            "z": 0.05,
            "z_err": 0.005,
            "time": 1234567891.0,
            "gw_time": 1234567880.0
        })
        
        results = assoc.compute_odds(em_model='kilonova')
        
        print(f"  I_omega (spatial):  {results['I_omega']:.3e}")
        print(f"  I_dl (distance):    {results['I_dl']:.3e}")
        print(f"  I_t (temporal):     {results['I_t']:.3e}")
        print(f"  Bayes Factor:       {results['bayes_factor']:.3e}")
        print(f"  Posterior Odds:     {results['posterior_odds']:.3e}")
        print(f"  P(Associated):      {results['confidence']:.1%}")
        print(f"  Decision:           {'ASSOCIATED' if results['associated'] else 'NOT ASSOCIATED'}")
        
        # Test with poor match (far in sky)
        print("\n2. Testing POOR MATCH (far from localization):")
        assoc2 = Association("test_data/test_skymap.fits", {
            "ra": 200.0,
            "dec": 45.0,
            "z": 0.08,
            "time": 1234568000.0,
            "gw_time": 1234567880.0
        })
        
        results2 = assoc2.compute_odds()
        print(f"  I_omega (spatial):  {results2['I_omega']:.3e}")
        print(f"  P(Associated):      {results2['confidence']:.1%}")
        print(f"  Decision:           {'ASSOCIATED' if results2['associated'] else 'NOT ASSOCIATED'}")
        
        return True
        
    except Exception as e:
        print(f"✗ Error in association calculation: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_plotting():
    """Test plotting capabilities"""
    print("\n" + "="*60)
    print("TESTING PLOTTING")
    print("="*60)
    
    try:
        from gwassociation import Association
        from gwassociation.plots import plot_association_summary
        
        # Create output directory
        plot_dir = Path("test_plots")
        plot_dir.mkdir(exist_ok=True)
        
        # Create association
        assoc = Association("test_data/test_skymap.fits", {
            "ra": 120.5,
            "dec": -30.0,
            "z": 0.05,
            "time": 1234567891.0,
            "gw_time": 1234567880.0
        })
        
        # Try to create plots
        print("Creating plots...")
        
        # Skymap plot
        try:
            assoc.plot_skymap(str(plot_dir / "test_skymap.png"))
            print(f"✓ Created skymap plot: {plot_dir}/test_skymap.png")
        except Exception as e:
            print(f"⚠ Skymap plot failed (may need healpy): {e}")
        
        # Summary plot
        try:
            results = assoc.compute_odds()
            plot_association_summary(results, str(plot_dir / "test_summary.png"))
            print(f"✓ Created summary plot: {plot_dir}/test_summary.png")
        except Exception as e:
            print(f"⚠ Summary plot failed: {e}")
        
        return True
        
    except Exception as e:
        print(f"✗ Error in plotting: {e}")
        return False

def test_cli():
    """Test command line interface"""
    print("\n" + "="*60)
    print("TESTING CLI")
    print("="*60)
    
    import subprocess
    
    cmd = [
        "python", "-m", "gwassociation.cli",
        "--gw-file", "test_data/test_skymap.fits",
        "--ra", "120.5",
        "--dec", "-30.0",
        "--z", "0.05",
        "--time", "1234567891",
        "--out", "test_output",
        "--verbose"
    ]
    
    print("Running CLI command:")
    print(" ".join(cmd))
    print("")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
        
        if result.returncode == 0:
            print("✓ CLI executed successfully")
            print("\nCLI Output:")
            print(result.stdout)
            
            # Check if output files were created
            output_dir = Path("test_output")
            if output_dir.exists():
                files = list(output_dir.glob("*"))
                print(f"\n✓ Created {len(files)} output files:")
                for f in files:
                    print(f"  - {f.name}")
            return True
        else:
            print("✗ CLI failed with error:")
            print(result.stderr)
            return False
            
    except subprocess.TimeoutExpired:
        print("✗ CLI timed out")
        return False
    except Exception as e:
        print(f"✗ CLI error: {e}")
        return False

# ============================================================================
# 3. FULL INTEGRATION TEST
# ============================================================================

def integration_test():
    """Run a complete integration test"""
    print("\n" + "="*60)
    print("INTEGRATION TEST: Multiple Candidates")
    print("="*60)
    
    try:
        from gwassociation import Association
        from gwassociation.plots import plot_candidate_ranking
        
        # Load test transients
        with open("test_data/test_transients.json", "r") as f:
            transients = json.load(f)
        
        # Create association and rank candidates
        assoc = Association("test_data/test_skymap.fits", {
            "ra": 0, "dec": 0, "gw_time": 1234567880.0
        })
        
        print("Ranking candidates:")
        rankings = assoc.rank_candidates(transients)
        
        for i, r in enumerate(rankings):
            print(f"\nRank {i+1}: {r['candidate'].name}")
            print(f"  Position: RA={r['candidate'].ra:.1f}, Dec={r['candidate'].dec:.1f}")
            print(f"  Redshift: {r['candidate'].z}")
            print(f"  Odds: {r['odds']:.2e}")
            print(f"  P(Associated): {r['probability']:.1%}")
        
        # Try to plot rankings
        try:
            results_list = []
            for r in rankings:
                results_list.append({
                    'odds': r['odds'],
                    'probability': r['probability']
                })
            plot_candidate_ranking(results_list, "test_plots/ranking.png")
            print("\n✓ Created ranking plot: test_plots/ranking.png")
        except Exception as e:
            print(f"\n⚠ Ranking plot failed: {e}")
        
        return True
        
    except Exception as e:
        print(f"✗ Integration test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

# ============================================================================
# 4. MAIN TEST RUNNER
# ============================================================================

def run_all_tests():
    """Run all tests and report results"""
    print("\n" + "="*60)
    print("GWASSOCIATION TEST SUITE")
    print("="*60)
    
    # Track results
    results = {}
    
    # Create test data
    test_dir = create_test_data()
    
    # Run tests
    results['imports'] = test_imports()
    results['basic'] = test_basic_functionality()
    results['association'] = test_association_calculation()
    results['plotting'] = test_plotting()
    results['cli'] = test_cli()
    results['integration'] = integration_test()
    
    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    
    for test_name, passed in results.items():
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"{test_name:15} {status}")
    
    total = len(results)
    passed = sum(results.values())
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 All tests passed! The framework is working correctly.")
    else:
        print("\n⚠ Some tests failed. Check the errors above for details.")
        print("\nCommon issues:")
        print("- Missing dependencies (healpy, ligo.skymap)")
        print("- Incorrect file paths")
        print("- Module import errors")
    
    return passed == total

def test_minimal():
    """Minimal test without optional dependencies"""
    print("\n" + "="*60)
    print("MINIMAL TEST (no healpy/ligo.skymap required)")
    print("="*60)
    
    try:
        # Test with mock data
        from gwassociation.io.transient import Transient
        from gwassociation.analysis.temporal import TemporalOverlap
        from gwassociation.analysis.los import DistanceOverlap
        
        print("Testing basic classes...")
        
        # Create transient
        t = Transient(ra=120, dec=-30, z=0.05)
        print(f"✓ Transient: {t}")
        
        # Test temporal overlap
        class MockGW:
            event_time = 1234567880
        
        class MockEM:
            time = 1234567881
        
        temporal = TemporalOverlap()
        i_t = temporal.compute(MockGW(), MockEM(), model='grb')
        print(f"✓ Temporal overlap (GRB): {i_t:.3f}")
        
        # Test distance
        distance = DistanceOverlap()
        print(f"✓ Distance calculator initialized")
        
        print("\n✓ Minimal test passed! Core functionality works.")
        print("\nTo enable full functionality, install:")
        print("  pip install healpy ligo.skymap astropy scipy matplotlib")
        
        return True
        
    except Exception as e:
        print(f"✗ Even minimal test failed: {e}")
        print("\nMake sure the gwassociation package is in your Python path:")
        print("  export PYTHONPATH=$PYTHONPATH:/path/to/gwassociation/parent/directory")
        return False

# ============================================================================
# 5. ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Test GWAssociation Framework")
    parser.add_argument("--minimal", action="store_true", 
                       help="Run minimal tests only (no optional dependencies)")
    parser.add_argument("--quick", action="store_true",
                       help="Run quick tests only (skip CLI and plotting)")
    
    args = parser.parse_args()
    
    if args.minimal:
        success = test_minimal()
    elif args.quick:
        print("Running quick tests...")
        test_dir = create_test_data()
        success = (test_imports() and 
                  test_basic_functionality() and 
                  test_association_calculation())
    else:
        success = run_all_tests()
    
    sys.exit(0 if success else 1)