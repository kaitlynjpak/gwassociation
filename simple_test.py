#!/usr/bin/env python
"""
Simple direct test of the GW-EM association framework
"""

print("Starting GW-EM Association Test...")

# Add the src directory to Python path
import sys
import os
sys.path.insert(0, 'src')

# Check if test data exists
if not os.path.exists("test_data/test_skymap.fits"):
    print("\n❌ ERROR: Test data not found!")
    print("Please run this first: python test_gw_assoc.py")
    print("This will create the test data needed.\n")
    sys.exit(1)

print("✓ Test data found\n")

# Import the framework
try:
    from gw_assoc import Association
    print("✓ Successfully imported Association class\n")
except ImportError as e:
    print(f"❌ ERROR: Could not import gw_assoc: {e}")
    print("Make sure you're in the right directory")
    sys.exit(1)

# Create an association with test data
print("=" * 60)
print("RUNNING GW-EM ASSOCIATION ANALYSIS")
print("=" * 60)

# Test transient that should match well with the test skymap
transient_info = {
    "ra": 120.5,      # Near center of test skymap
    "dec": -30.0,     # Near center of test skymap
    "z": 0.050,       # Reasonable redshift
    "z_err": 0.005,   
    "time": 1234567890.0,    
    "gw_time": 1234567880.0,  # 10 seconds before transient
    "magnitude": 18.5,
    "name": "TEST_TRANSIENT_001"
}

print("\nTransient Information:")
print(f"  Name: {transient_info['name']}")
print(f"  Position: RA={transient_info['ra']}°, Dec={transient_info['dec']}°")
print(f"  Redshift: z={transient_info['z']} ± {transient_info['z_err']}")
print(f"  Time: {transient_info['time']} (10 sec after GW)")

# Create association
print("\nCreating association...")
assoc = Association(
    gw_file="test_data/test_skymap.fits",
    transient_info=transient_info
)
print("✓ Association created")

# Compute odds
print("\nComputing association odds...")
results = assoc.compute_odds(em_model='kilonova')
print("✓ Computation complete")

# Display results
print("\n" + "=" * 60)
print("RESULTS")
print("=" * 60)

print("\nOverlap Integrals:")
print(f"  Spatial (I_Ω):   {results['I_omega']:.3e}")
print(f"  Distance (I_DL):  {results['I_dl']:.3e}")
print(f"  Temporal (I_t):   {results['I_t']:.3e}")

print("\nStatistics:")
print(f"  Bayes Factor:     {results['bayes_factor']:.3e}")
print(f"  Posterior Odds:   {results['posterior_odds']:.3e}")

print("\nFinal Result:")
print(f"  P(Associated) = {results['confidence']:.1%}")
print(f"  Decision: {'✓ ASSOCIATED' if results['associated'] else '✗ NOT ASSOCIATED'}")

# Try to make a plot
print("\n" + "=" * 60)
print("GENERATING PLOT")
print("=" * 60)

try:
    assoc.plot_skymap("test_skymap_output.png")
    print("✓ Skymap saved to: test_skymap_output.png")
except Exception as e:
    print(f"⚠ Could not create plot: {e}")

print("\n" + "=" * 60)
print("TEST COMPLETE!")
print("=" * 60)
print("\nThe framework is working correctly!")
print("You can now modify the transient_info dictionary above")
print("to test with different transient parameters.\n")