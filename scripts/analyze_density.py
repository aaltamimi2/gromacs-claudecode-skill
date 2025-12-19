#!/usr/bin/env python3
"""
Density analysis for GROMACS trajectories

Calculates system density at regular intervals from MD trajectories.
Uses MDAnalysis for trajectory processing and provides detailed statistics.

Usage:
    python analyze_density.py <topology.tpr> <trajectory.xtc> [options]

Example:
    python analyze_density.py npt.tpr npt.xtc -o density.dat -i 2.0 --plot
"""

import sys
import argparse
from pathlib import Path
import numpy as np

# Check for MDAnalysis
try:
    import MDAnalysis as mda
    from MDAnalysis.lib.log import ProgressBar
except ImportError:
    print("ERROR: MDAnalysis not found. Install with:")
    print("  pip install MDAnalysis")
    print("  or: conda install -c conda-forge MDAnalysis")
    sys.exit(1)


def calculate_density(topology, trajectory, output_file="density_vs_time.dat",
                      interval_ns=2.0, skip_frames=0, plot=False):
    """
    Calculate density of the system at regular intervals.

    Parameters:
    -----------
    topology : str
        Path to topology file (.tpr, .gro, .pdb)
    trajectory : str
        Path to trajectory file (.xtc, .trr)
    output_file : str
        Output file for density data
    interval_ns : float
        Time interval in nanoseconds for density calculation
    skip_frames : int
        Number of initial frames to skip (for equilibration)
    plot : bool
        Generate matplotlib plot if True

    Returns:
    --------
    times : list
        Time points in nanoseconds
    densities : list
        Densities in g/cm³
    """

    print("=" * 60)
    print("GROMACS Density Analysis")
    print("=" * 60)
    print(f"Topology:    {topology}")
    print(f"Trajectory:  {trajectory}")
    print(f"Output:      {output_file}")
    print(f"Interval:    {interval_ns} ns")
    if skip_frames > 0:
        print(f"Skipping:    {skip_frames} initial frames")
    print("=" * 60)

    # Load trajectory
    try:
        u = mda.Universe(topology, trajectory)
    except Exception as e:
        print(f"\nERROR: Failed to load trajectory")
        print(f"  {e}")
        sys.exit(1)

    # Check trajectory is not empty
    if len(u.trajectory) == 0:
        print("\nERROR: Trajectory is empty")
        sys.exit(1)

    # Get trajectory parameters
    dt_ps = u.trajectory.dt  # Time step in ps
    interval_ps = interval_ns * 1000  # Convert ns to ps
    frame_interval = max(1, int(interval_ps / dt_ps))

    total_frames = len(u.trajectory)
    total_time_ns = u.trajectory.totaltime / 1000

    print(f"\nTrajectory Information:")
    print(f"  Time step:     {dt_ps:.3f} ps")
    print(f"  Total frames:  {total_frames}")
    print(f"  Total time:    {total_time_ns:.2f} ns")
    print(f"  Frame interval: {frame_interval} frames ({interval_ns} ns)")

    if skip_frames >= total_frames:
        print(f"\nERROR: skip_frames ({skip_frames}) >= total frames ({total_frames})")
        sys.exit(1)

    # Storage for results
    times = []
    densities = []
    volumes = []

    # Calculate density at each interval
    print(f"\nCalculating density...")
    print(f"{'Time (ns)':>10} {'Density (g/cm³)':>18} {'Volume (nm³)':>15}")
    print("-" * 45)

    try:
        for ts in ProgressBar(u.trajectory[skip_frames::frame_interval],
                             verbose=True, desc="Processing"):
            time_ns = ts.time / 1000  # Convert ps to ns

            # Get box volume in Å³
            box_volume_angstrom3 = ts.volume

            # Convert to nm³ for display
            box_volume_nm3 = box_volume_angstrom3 / 1000

            # Calculate total mass in atomic mass units (amu)
            total_mass_amu = u.atoms.total_mass()

            # Calculate density in g/cm³
            # Formula: ρ = (M / V) * (N_A / N_A) * (1 g/mol / 1 amu) * (1 Å³ / 1e-24 cm³)
            # Simplified: ρ = M_amu / V_Å³ * 1.66054
            # where 1.66054 = 1/(N_A * 1e-24) and converts amu/Å³ to g/cm³
            density_g_cm3 = total_mass_amu / box_volume_angstrom3 * 1.66054

            times.append(time_ns)
            densities.append(density_g_cm3)
            volumes.append(box_volume_nm3)

            print(f"{time_ns:10.2f} {density_g_cm3:18.6f} {box_volume_nm3:15.3f}")

    except KeyboardInterrupt:
        print("\n\nAnalysis interrupted by user")
        if len(densities) == 0:
            print("No data collected. Exiting.")
            sys.exit(1)
        print(f"Collected {len(densities)} data points. Continuing with analysis...")

    except Exception as e:
        print(f"\nERROR during trajectory analysis:")
        print(f"  {e}")
        sys.exit(1)

    # Convert to numpy arrays for statistics
    times = np.array(times)
    densities = np.array(densities)
    volumes = np.array(volumes)

    # Save results
    print(f"\n{'='*60}")
    print(f"Saving results to {output_file}")
    try:
        with open(output_file, 'w') as f:
            f.write("# Density analysis from GROMACS trajectory\n")
            f.write(f"# Topology: {topology}\n")
            f.write(f"# Trajectory: {trajectory}\n")
            f.write(f"# Interval: {interval_ns} ns\n")
            f.write(f"# Skipped frames: {skip_frames}\n")
            f.write("#\n")
            f.write("# Time (ns)    Density (g/cm³)    Volume (nm³)\n")
            for t, d, v in zip(times, densities, volumes):
                f.write(f"{t:12.4f}    {d:16.8f}    {v:14.6f}\n")
        print(f"✓ Data saved successfully")
    except Exception as e:
        print(f"ERROR: Failed to save output file")
        print(f"  {e}")
        sys.exit(1)

    # Calculate statistics
    mean_density = np.mean(densities)
    std_density = np.std(densities)
    sem_density = std_density / np.sqrt(len(densities))  # Standard error

    mean_volume = np.mean(volumes)
    std_volume = np.std(volumes)

    print(f"\n{'='*60}")
    print("DENSITY STATISTICS")
    print("=" * 60)
    print(f"  Mean:              {mean_density:.6f} ± {std_density:.6f} g/cm³")
    print(f"  Std Error (SEM):   {sem_density:.6f} g/cm³")
    print(f"  Min:               {np.min(densities):.6f} g/cm³")
    print(f"  Max:               {np.max(densities):.6f} g/cm³")
    print(f"  Data points:       {len(densities)}")
    print(f"\nVOLUME STATISTICS")
    print("-" * 60)
    print(f"  Mean:              {mean_volume:.4f} ± {std_volume:.4f} nm³")
    print(f"  Min:               {np.min(volumes):.4f} nm³")
    print(f"  Max:               {np.max(volumes):.4f} nm³")

    # Check for drift (compare first and last 20%)
    n_points = len(densities)
    n_check = max(1, n_points // 5)  # 20% of points

    first_mean = np.mean(densities[:n_check])
    last_mean = np.mean(densities[-n_check:])
    drift = last_mean - first_mean
    drift_percent = (drift / first_mean) * 100

    print(f"\nDRIFT ANALYSIS")
    print("-" * 60)
    print(f"  First 20% mean:    {first_mean:.6f} g/cm³")
    print(f"  Last 20% mean:     {last_mean:.6f} g/cm³")
    print(f"  Drift:             {drift:+.6f} g/cm³ ({drift_percent:+.2f}%)")

    if abs(drift_percent) > 1.0:
        print(f"  ⚠️  WARNING: Density drift > 1%, system may not be equilibrated")
    else:
        print(f"  ✓ Density stable (drift < 1%)")

    print("=" * 60)

    # Generate plot if requested
    if plot:
        try:
            import matplotlib.pyplot as plt

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

            # Density plot
            ax1.plot(times, densities, 'b-', linewidth=1.5, label='Density')
            ax1.axhline(mean_density, color='r', linestyle='--',
                       label=f'Mean: {mean_density:.4f} g/cm³')
            ax1.fill_between(times,
                            mean_density - std_density,
                            mean_density + std_density,
                            alpha=0.3, color='r', label='±1 σ')
            ax1.set_ylabel('Density (g/cm³)', fontsize=12)
            ax1.legend(loc='best')
            ax1.grid(True, alpha=0.3)

            # Volume plot
            ax2.plot(times, volumes, 'g-', linewidth=1.5, label='Volume')
            ax2.axhline(mean_volume, color='orange', linestyle='--',
                       label=f'Mean: {mean_volume:.2f} nm³')
            ax2.set_xlabel('Time (ns)', fontsize=12)
            ax2.set_ylabel('Volume (nm³)', fontsize=12)
            ax2.legend(loc='best')
            ax2.grid(True, alpha=0.3)

            plt.suptitle(f'Density Analysis: {Path(trajectory).stem}', fontsize=14)
            plt.tight_layout()

            plot_file = output_file.replace('.dat', '.png')
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            print(f"\n✓ Plot saved to {plot_file}")
            plt.close()

        except ImportError:
            print("\n⚠️  matplotlib not found. Install for plotting:")
            print("  pip install matplotlib")
        except Exception as e:
            print(f"\n⚠️  Failed to generate plot: {e}")

    return times, densities, volumes


def main():
    """Main entry point with argument parsing."""
    parser = argparse.ArgumentParser(
        description="Analyze density from GROMACS trajectory",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (2 ns intervals)
  python analyze_density.py npt.tpr npt.xtc

  # Custom interval and output
  python analyze_density.py npt.tpr npt.xtc -o my_density.dat -i 1.0

  # Skip first 1000 frames and plot
  python analyze_density.py npt.tpr npt.xtc --skip 1000 --plot

  # High-resolution sampling (0.5 ns intervals)
  python analyze_density.py npt.tpr npt.xtc -i 0.5 --plot
        """
    )

    parser.add_argument('topology',
                       help='Topology file (.tpr, .gro, .pdb)')
    parser.add_argument('trajectory',
                       help='Trajectory file (.xtc, .trr)')
    parser.add_argument('-o', '--output',
                       default='density_vs_time.dat',
                       help='Output file (default: density_vs_time.dat)')
    parser.add_argument('-i', '--interval',
                       type=float,
                       default=2.0,
                       help='Sampling interval in ns (default: 2.0)')
    parser.add_argument('--skip',
                       type=int,
                       default=0,
                       help='Number of initial frames to skip (default: 0)')
    parser.add_argument('--plot',
                       action='store_true',
                       help='Generate density vs time plot (requires matplotlib)')

    args = parser.parse_args()

    # Validate inputs
    if not Path(args.topology).exists():
        print(f"ERROR: Topology file not found: {args.topology}")
        sys.exit(1)

    if not Path(args.trajectory).exists():
        print(f"ERROR: Trajectory file not found: {args.trajectory}")
        sys.exit(1)

    if args.interval <= 0:
        print(f"ERROR: Interval must be positive, got {args.interval}")
        sys.exit(1)

    if args.skip < 0:
        print(f"ERROR: Skip frames must be non-negative, got {args.skip}")
        sys.exit(1)

    # Run analysis
    try:
        times, densities, volumes = calculate_density(
            args.topology,
            args.trajectory,
            args.output,
            args.interval,
            args.skip,
            args.plot
        )
    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    print("\n✓ Analysis complete!\n")


if __name__ == "__main__":
    main()
