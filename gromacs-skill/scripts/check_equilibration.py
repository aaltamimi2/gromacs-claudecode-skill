#!/usr/bin/env python3
"""
check_equilibration.py - Validate GROMACS equilibration runs

Usage:
    python check_equilibration.py nvt.edr npt.edr [--target-temp 300] [--target-density 1000]
"""

import subprocess
import argparse
import sys
from pathlib import Path


def extract_energy_data(edr_file, term):
    """Extract data from .edr using gmx energy."""
    cmd = f"echo '{term}' | gmx energy -f {edr_file} -o /dev/stdout 2>/dev/null"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    lines = [l for l in result.stdout.strip().split('\n') 
             if l and not l.startswith('#') and not l.startswith('@')]
    
    if not lines:
        return None, None, None
    
    values = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 2:
            try:
                values.append(float(parts[1]))
            except ValueError:
                continue
    
    if not values:
        return None, None, None
    
    mean = sum(values) / len(values)
    variance = sum((x - mean) ** 2 for x in values) / len(values)
    std = variance ** 0.5
    
    # Check for drift (compare first and last 10%)
    n = len(values)
    first_mean = sum(values[:max(1, n//10)]) / max(1, n//10)
    last_mean = sum(values[-max(1, n//10):]) / max(1, n//10)
    drift = last_mean - first_mean
    
    return mean, std, drift


def check_nvt(edr_file, target_temp=300, tolerance=5):
    """Check NVT equilibration."""
    print(f"\n=== NVT Check: {edr_file} ===")
    
    mean, std, drift = extract_energy_data(edr_file, "Temperature")
    
    if mean is None:
        print("  ERROR: Could not extract temperature data")
        return False
    
    print(f"  Temperature: {mean:.1f} ± {std:.1f} K")
    print(f"  Drift: {drift:+.2f} K")
    
    passed = True
    
    if abs(mean - target_temp) > tolerance:
        print(f"  WARNING: Temperature deviates from target ({target_temp} K)")
        passed = False
    
    if abs(drift) > tolerance:
        print(f"  WARNING: Temperature drift detected")
        passed = False
    
    if passed:
        print("  ✓ NVT equilibration looks good")
    
    return passed


def check_npt(edr_file, target_temp=300, target_density=1000, temp_tol=5, dens_tol=50):
    """Check NPT equilibration."""
    print(f"\n=== NPT Check: {edr_file} ===")
    
    passed = True
    
    # Temperature
    mean_t, std_t, drift_t = extract_energy_data(edr_file, "Temperature")
    if mean_t is not None:
        print(f"  Temperature: {mean_t:.1f} ± {std_t:.1f} K (drift: {drift_t:+.2f} K)")
        if abs(mean_t - target_temp) > temp_tol:
            print(f"  WARNING: Temperature deviates from target ({target_temp} K)")
            passed = False
    
    # Density
    mean_d, std_d, drift_d = extract_energy_data(edr_file, "Density")
    if mean_d is not None:
        print(f"  Density: {mean_d:.1f} ± {std_d:.1f} kg/m³ (drift: {drift_d:+.2f})")
        if abs(mean_d - target_density) > dens_tol:
            print(f"  WARNING: Density deviates from target ({target_density} kg/m³)")
            passed = False
    
    # Pressure (for reference - fluctuations are normal)
    mean_p, std_p, drift_p = extract_energy_data(edr_file, "Pressure")
    if mean_p is not None:
        print(f"  Pressure: {mean_p:.1f} ± {std_p:.1f} bar (large fluctuations normal)")
    
    # Potential energy
    mean_e, std_e, drift_e = extract_energy_data(edr_file, "Potential")
    if mean_e is not None:
        rel_drift = abs(drift_e / mean_e) * 100 if mean_e != 0 else 0
        print(f"  Potential energy: {mean_e:.0f} ± {std_e:.0f} kJ/mol")
        if rel_drift > 1:
            print(f"  WARNING: Potential energy drift ({rel_drift:.1f}%)")
    
    if passed:
        print("  ✓ NPT equilibration looks good")
    
    return passed


def main():
    parser = argparse.ArgumentParser(description="Check GROMACS equilibration")
    parser.add_argument("edr_files", nargs="+", help="Energy files to check")
    parser.add_argument("--target-temp", type=float, default=300, help="Target temperature (K)")
    parser.add_argument("--target-density", type=float, default=1000, help="Target density (kg/m³)")
    args = parser.parse_args()
    
    all_passed = True
    
    for edr in args.edr_files:
        if not Path(edr).exists():
            print(f"ERROR: File not found: {edr}")
            all_passed = False
            continue
        
        name = Path(edr).stem.lower()
        
        if 'nvt' in name:
            passed = check_nvt(edr, args.target_temp)
        elif 'npt' in name or 'md' in name:
            passed = check_npt(edr, args.target_temp, args.target_density)
        else:
            # Default to NPT check
            passed = check_npt(edr, args.target_temp, args.target_density)
        
        all_passed = all_passed and passed
    
    print("\n" + "=" * 40)
    if all_passed:
        print("Overall: ✓ All checks passed")
        sys.exit(0)
    else:
        print("Overall: ✗ Some checks failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
