#!/usr/bin/env python3
"""
Conformational sampling from MD trajectories using Rg-SASA grid

Samples conformers across radius of gyration (Rg) and solvent-accessible
surface area (SASA) space for subsequent DFT calculations (e.g., COSMO-RS).

Usage:
    python conformational_sampling.py npt.xtc npt.gro -g 10 --plot

Example:
    python conformational_sampling.py npt.xtc npt.gro --grid-size 10 \
        --threshold 0.03 --output conformers.dat --plot
"""

import sys
import argparse
from pathlib import Path
import numpy as np

# Check for required packages
try:
    import mdtraj as md
except ImportError:
    print("ERROR: MDTraj not found. Install with:")
    print("  pip install mdtraj")
    print("  or: conda install -c conda-forge mdtraj")
    sys.exit(1)

try:
    import matplotlib
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


def normalize(vector):
    """Normalize vector to [0, 1] range."""
    v_min = vector.min()
    v_max = vector.max()
    if v_max - v_min == 0:
        return np.zeros_like(vector)
    return (vector - v_min) / (v_max - v_min)


def generate_grid_coordinates(n_points):
    """
    Generate n×n grid of coordinates in [0, 1]×[0, 1] space.

    Parameters:
    -----------
    n_points : int
        Number of grid points along each dimension

    Returns:
    --------
    coordinates : ndarray
        Array of (x, y) grid coordinates, shape (n*n, 2)
    """
    coordinates = np.zeros((n_points * n_points, 2))
    for i in range(n_points):
        for j in range(n_points):
            coordinates[i * n_points + j] = (i / (n_points - 1), j / (n_points - 1))
    return coordinates


def find_nearest_conformers(grid_coords, rg_norm, sasa_norm, threshold=0.03):
    """
    Find nearest conformer to each grid point within threshold distance.

    Parameters:
    -----------
    grid_coords : ndarray
        Grid coordinates, shape (n_grid, 2)
    rg_norm : ndarray
        Normalized Rg values for all frames
    sasa_norm : ndarray
        Normalized SASA values for all frames
    threshold : float
        Maximum distance for a grid point to be considered "filled"

    Returns:
    --------
    selected_indices : list
        Frame indices of selected conformers
    grid_points : list
        Grid coordinates that have conformers
    conformer_coords : list
        Actual (Rg, SASA) coordinates of selected conformers
    """
    selected_indices = []
    grid_points = []
    conformer_coords = []

    for x, y in grid_coords:
        # Calculate distance from grid point to all conformers
        distances = np.sqrt((rg_norm - x)**2 + (sasa_norm - y)**2)
        nearest_idx = np.argmin(distances)

        # Only include if within threshold
        if distances[nearest_idx] <= threshold:
            # Avoid duplicates
            if nearest_idx not in selected_indices:
                selected_indices.append(nearest_idx)
                grid_points.append((x, y))
                conformer_coords.append((rg_norm[nearest_idx], sasa_norm[nearest_idx]))

    return selected_indices, grid_points, conformer_coords


def analyze_conformers(trajectory, topology, grid_size=10, threshold=0.03,
                       output_file="conformers.dat", plot=False):
    """
    Perform conformational sampling from MD trajectory.

    Parameters:
    -----------
    trajectory : str
        Path to trajectory file (.xtc, .trr)
    topology : str
        Path to topology file (.gro, .pdb, .tpr)
    grid_size : int
        Number of grid points along each dimension (total = grid_size²)
    threshold : float
        Maximum distance in normalized space for grid point matching
    output_file : str
        Output file for selected frame indices
    plot : bool
        Generate matplotlib plots

    Returns:
    --------
    selected_frames : list
        Frame indices to extract for DFT
    """

    print("=" * 70)
    print("CONFORMATIONAL SAMPLING - Rg-SASA Grid")
    print("=" * 70)
    print(f"Trajectory:  {trajectory}")
    print(f"Topology:    {topology}")
    print(f"Grid size:   {grid_size}×{grid_size} = {grid_size**2} points")
    print(f"Threshold:   {threshold} (normalized units)")
    print("\n⚠️  Reminder: Use a polymer-only trajectory for grid/umbrella sampling.")
    print('   Example: echo "Polymer" | gmx trjconv -f npt.xtc -s npt.gro -o polymer_only.xtc')
    print("=" * 70)

    # Load trajectory
    print("\nLoading trajectory...")
    try:
        traj = md.load(trajectory, top=topology)
    except Exception as e:
        print(f"ERROR: Failed to load trajectory")
        print(f"  {e}")
        sys.exit(1)

    n_frames = traj.n_frames
    print(f"✓ Loaded {n_frames} frames")

    # Calculate Rg and SASA
    print("\nCalculating Rg and SASA...")
    try:
        sasa = md.shrake_rupley(traj).sum(axis=1)  # Total SASA per frame
        rg = md.compute_rg(traj)
    except Exception as e:
        print(f"ERROR: Failed to calculate properties")
        print(f"  {e}")
        sys.exit(1)

    print(f"✓ Rg range:   {rg.min():.3f} - {rg.max():.3f} nm")
    print(f"✓ SASA range: {sasa.min():.1f} - {sasa.max():.1f} nm²")

    # Normalize
    rg_norm = normalize(rg)
    sasa_norm = normalize(sasa)

    # Generate grid
    print(f"\nGenerating {grid_size}×{grid_size} grid...")
    grid_coords = generate_grid_coordinates(grid_size)

    # Find nearest conformers
    print(f"Finding conformers within threshold {threshold}...")
    selected_indices, grid_points, conformer_coords = find_nearest_conformers(
        grid_coords, rg_norm, sasa_norm, threshold
    )

    n_selected = len(selected_indices)
    coverage = (n_selected / len(grid_coords)) * 100

    print(f"\n✓ Selected {n_selected} conformers")
    print(f"  Grid coverage: {coverage:.1f}%")

    # Find extreme conformers
    lowest_rg_indices = np.argsort(rg)[:5]
    lowest_sasa_indices = np.argsort(sasa)[:5]
    highest_rg_indices = np.argsort(rg)[-5:][::-1]
    highest_sasa_indices = np.argsort(sasa)[-5:][::-1]

    # Combine all selections (remove duplicates)
    all_selected = set(selected_indices)
    all_selected.update(lowest_rg_indices)
    all_selected.update(lowest_sasa_indices)
    all_selected.update(highest_rg_indices)
    all_selected.update(highest_sasa_indices)

    final_indices = sorted(all_selected)

    print(f"\n✓ Total conformers (grid + extremes): {len(final_indices)}")

    # Capture frame times (ps) for extraction script
    frame_times_ps = traj.time

    # Save results (original format)
    print(f"\nSaving results to {output_file}...")
    try:
        with open(output_file, 'w') as f:
            f.write("# Conformational sampling results\n")
            f.write(f"# Trajectory: {trajectory}\n")
            f.write(f"# Grid size: {grid_size}×{grid_size}\n")
            f.write(f"# Threshold: {threshold}\n")
            f.write(f"# Total selected: {len(final_indices)}\n")
            f.write("#\n")
            f.write("# Frame    Rg(nm)    SASA(nm²)    Rg_norm    SASA_norm    Source\n")

            for idx in final_indices:
                source = []
                if idx in selected_indices:
                    source.append("grid")
                if idx in lowest_rg_indices:
                    source.append("low_Rg")
                if idx in lowest_sasa_indices:
                    source.append("low_SASA")
                if idx in highest_rg_indices:
                    source.append("high_Rg")
                if idx in highest_sasa_indices:
                    source.append("high_SASA")

                source_str = "+".join(source)

                f.write(f"{idx:8d}  {rg[idx]:8.4f}  {sasa[idx]:10.2f}  "
                       f"{rg_norm[idx]:9.4f}  {sasa_norm[idx]:10.4f}  {source_str}\n")

        print(f"✓ Saved to {output_file}")
    except Exception as e:
        print(f"ERROR: Failed to save output")
        print(f"  {e}")
        sys.exit(1)

    # Save CSV metadata for easy cross-referencing with extracted PDB files
    csv_file = output_file.replace('.dat', '_metadata.csv')
    print(f"Saving metadata CSV to {csv_file}...")
    try:
        import csv
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['conformer_id', 'conformer_file', 'frame', 'time_ps',
                           'Rg_nm', 'SASA_nm2', 'Rg_norm', 'SASA_norm', 'source'])

            for i, idx in enumerate(final_indices):
                conformer_id = i + 1
                conformer_file = f"conformer_{conformer_id:03d}_frame{idx}.pdb"
                time_ps = float(frame_times_ps[idx])

                source = []
                if idx in selected_indices:
                    source.append("grid")
                if idx in lowest_rg_indices:
                    source.append("low_Rg")
                if idx in lowest_sasa_indices:
                    source.append("low_SASA")
                if idx in highest_rg_indices:
                    source.append("high_Rg")
                if idx in highest_sasa_indices:
                    source.append("high_SASA")
                source_str = "+".join(source)

                writer.writerow([conformer_id, conformer_file, idx, f"{time_ps:.3f}",
                               f"{rg[idx]:.4f}", f"{sasa[idx]:.2f}",
                               f"{rg_norm[idx]:.4f}", f"{sasa_norm[idx]:.4f}",
                               source_str])

        print(f"✓ Saved metadata to {csv_file}")
        print(f"  Use this file to cross-reference conformer PDB files with Rg-SASA values")
    except Exception as e:
        print(f"WARNING: Failed to save CSV metadata: {e}")

    # Generate trjconv extraction script
    extract_script = output_file.replace('.dat', '_extract.sh')
    print(f"\nGenerating extraction script: {extract_script}")

    try:
        with open(extract_script, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# Extract selected conformers using gmx trjconv\n")
            f.write(f"# Generated from: {trajectory}\n")
            f.write(f"# IMPORTANT: Only polymer atoms should be extracted (not solvent!)\n\n")

            f.write("TRAJ=\"" + trajectory + "\"\n")
            f.write("TOP=\"" + topology + "\"\n")
            f.write("OUTDIR=\"conformers\"\n\n")

            f.write("mkdir -p ${OUTDIR}\n\n")

            f.write("# Detect Polymer group automatically (override with POLYMER_GROUP)\n")
            f.write("POLYMER_GROUP=${POLYMER_GROUP:-}\n")
            f.write("detect_polymer_group() {\n")
            f.write("    local top_file=\"$1\"\n")
            f.write("    if ! command -v gmx >/dev/null 2>&1; then\n")
            f.write("        return\n")
            f.write("    fi\n")
            f.write("    local ndx_output\n")
            f.write("    ndx_output=$(gmx make_ndx -f \"$top_file\" -quiet << EOF\nq\nEOF\n)\n")
            f.write("    echo \"$ndx_output\" | awk 'toupper($0) ~ /POLYMER/ && $1 ~ /^[0-9]+$/ {print $1; exit}'\n")
            f.write("}\n\n")

            f.write('if [ -z "$POLYMER_GROUP" ]; then\n')
            f.write('    POLYMER_GROUP=$(detect_polymer_group "${TOP}")\n')
            f.write('fi\n')
            f.write('if [ -z "$POLYMER_GROUP" ]; then\n')
            f.write('    echo "⚠️  Could not auto-detect Polymer group; defaulting to group 1."\n')
            f.write('    POLYMER_GROUP=1\n')
            f.write('else\n')
            f.write('    echo "Using Polymer group index: ${POLYMER_GROUP}"\n')
            f.write('fi\n\n')

            f.write("echo \"Extracting conformers using group ${POLYMER_GROUP} (polymer only)...\"\n")
            f.write("echo \"\"\n\n")

            for i, idx in enumerate(final_indices):
                time_ps = float(frame_times_ps[idx])
                f.write(f"# Conformer {i+1}: Frame {idx} (t = {time_ps:.3f} ps, Rg = {rg[idx]:.3f} nm, SASA = {sasa[idx]:.1f} nm²)\n")
                f.write(f"echo ${{POLYMER_GROUP}} | gmx trjconv -f ${{TRAJ}} -s ${{TOP}} "
                       f"-dump {time_ps:.3f} -o ${{OUTDIR}}/conformer_{i+1:03d}_frame{idx}.pdb > /dev/null 2>&1\n\n")

            f.write(f"echo \"\"\n")
            f.write(f"echo \"✓ Extracted {len(final_indices)} conformers to ${{OUTDIR}}/\"\n")
            f.write(f"echo \"✓ PDB files contain POLYMER ONLY (no solvent)\"\n")
            f.write(f"echo \"✓ Metadata saved to conformers_metadata.csv\"\n")
            f.write(f"echo \"\"\n")
            f.write(f"echo \"Next steps:\"\n")
            f.write(f"echo \"  1. Review conformers_metadata.csv to see Rg-SASA values for each conformer\"\n")
            f.write(f"echo \"  2. Open PDB files in GaussView to create .gjf files\"\n")
            f.write(f"echo \"  3. Run prepare_gaussian.py to add COSMO-RS commands\"\n")

        Path(extract_script).chmod(0o755)
        print(f"✓ Extraction script saved")
        print(f"\nTo extract conformers, run:")
        print(f"  ./{extract_script}")
        print(f"\nOr specify polymer group:")
        print(f"  POLYMER_GROUP=1 ./{extract_script}")

    except Exception as e:
        print(f"WARNING: Failed to create extraction script: {e}")

    # Print statistics
    print("\n" + "=" * 70)
    print("SELECTED CONFORMERS")
    print("=" * 70)

    print(f"\n=== Grid-sampled conformers: {len(selected_indices)} ===")
    for idx in selected_indices[:10]:  # Show first 10
        print(f"  Frame {idx:6d}: Rg = {rg[idx]:6.3f} nm, SASA = {sasa[idx]:8.1f} nm²")
    if len(selected_indices) > 10:
        print(f"  ... and {len(selected_indices) - 10} more")

    print(f"\n=== 5 Lowest Rg conformers ===")
    for idx in lowest_rg_indices:
        print(f"  Frame {idx:6d}: Rg = {rg[idx]:6.3f} nm, SASA = {sasa[idx]:8.1f} nm²")

    print(f"\n=== 5 Lowest SASA conformers ===")
    for idx in lowest_sasa_indices:
        print(f"  Frame {idx:6d}: SASA = {sasa[idx]:8.1f} nm², Rg = {rg[idx]:6.3f} nm")

    print("=" * 70)

    # Generate plots if requested
    if plot:
        if not MATPLOTLIB_AVAILABLE:
            print("\n⚠️  matplotlib not available. Install to enable plotting:")
            print("  pip install matplotlib")
        else:
            print("\nGenerating plots...")
            plot_results(rg, sasa, rg_norm, sasa_norm, selected_indices,
                        grid_points, conformer_coords, lowest_rg_indices,
                        lowest_sasa_indices, highest_rg_indices, highest_sasa_indices,
                        final_indices, output_file)

    return final_indices


def plot_results(rg, sasa, rg_norm, sasa_norm, selected_indices,
                grid_points, conformer_coords, lowest_rg, lowest_sasa,
                highest_rg, highest_sasa, final_indices, output_file):
    """Generate analysis plots with conformer annotations."""

    matplotlib.rcParams['font.family'] = 'Arial'
    fig = plt.figure(figsize=(15, 10))

    # 1. Rg histogram
    ax1 = plt.subplot(2, 3, 1)
    ax1.hist(rg, bins=40, edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Radius of Gyration (nm)', fontsize=11)
    ax1.set_ylabel('Count', fontsize=11)
    ax1.set_title('Rg Distribution', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)

    # 2. SASA histogram
    ax2 = plt.subplot(2, 3, 2)
    ax2.hist(sasa, bins=40, edgecolor='black', alpha=0.7)
    ax2.set_xlabel('SASA (nm²)', fontsize=11)
    ax2.set_ylabel('Count', fontsize=11)
    ax2.set_title('SASA Distribution', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    # 3. Normalized scatter
    ax3 = plt.subplot(2, 3, 3)
    ax3.scatter(rg_norm, sasa_norm, alpha=0.3, s=10, label='All frames')
    ax3.set_xlabel('Normalized Rg', fontsize=11)
    ax3.set_ylabel('Normalized SASA', fontsize=11)
    ax3.set_title('Normalized Rg-SASA Space', fontsize=12, fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # 4. Grid sampling visualization with conformer numbers
    ax4 = plt.subplot(2, 3, 4)
    # All frames
    ax4.scatter(rg_norm, sasa_norm, c='lightblue', alpha=0.3, s=10, label='All frames')
    # Grid points
    if grid_points:
        grid_x, grid_y = zip(*grid_points)
        ax4.scatter(grid_x, grid_y, c='red', s=30, marker='s', alpha=0.6, label='Grid points')
    # Selected conformers with annotations
    if conformer_coords:
        conf_x, conf_y = zip(*conformer_coords)
        ax4.scatter(conf_x, conf_y, c='green', s=50, marker='x', linewidths=2,
                   label='Selected conformers')

    # Annotate all final conformers with their IDs
    for i, idx in enumerate(final_indices):
        conformer_id = i + 1
        # Only annotate a subset to avoid overcrowding (every 3rd conformer for large sets)
        if len(final_indices) <= 30 or i % 3 == 0:
            ax4.annotate(f'{conformer_id}',
                        (rg_norm[idx], sasa_norm[idx]),
                        xytext=(3, 3), textcoords='offset points',
                        fontsize=6, alpha=0.7, color='darkgreen')

    ax4.set_xlabel('Normalized Rg', fontsize=11)
    ax4.set_ylabel('Normalized SASA', fontsize=11)
    ax4.set_title('Grid Sampling (with conformer IDs)', fontsize=12, fontweight='bold')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    # 5. Extreme conformers with annotations
    ax5 = plt.subplot(2, 3, 5)
    ax5.scatter(rg_norm, sasa_norm, c='lightgray', alpha=0.3, s=10, label='All frames')
    # Lowest Rg
    ax5.scatter(rg_norm[lowest_rg], sasa_norm[lowest_rg], c='blue', s=100,
               marker='o', edgecolors='black', linewidths=1.5, label='Lowest Rg')
    # Lowest SASA
    ax5.scatter(rg_norm[lowest_sasa], sasa_norm[lowest_sasa], c='orange', s=100,
               marker='s', edgecolors='black', linewidths=1.5, label='Lowest SASA')
    # Highest Rg
    ax5.scatter(rg_norm[highest_rg], sasa_norm[highest_rg], c='purple', s=100,
               marker='^', edgecolors='black', linewidths=1.5, label='Highest Rg')
    # Highest SASA
    ax5.scatter(rg_norm[highest_sasa], sasa_norm[highest_sasa], c='red', s=100,
               marker='v', edgecolors='black', linewidths=1.5, label='Highest SASA')

    # Annotate extreme conformers with their IDs
    for idx in list(lowest_rg) + list(lowest_sasa) + list(highest_rg) + list(highest_sasa):
        if idx in final_indices:
            conformer_id = final_indices.index(idx) + 1
            ax5.annotate(f'{conformer_id}',
                        (rg_norm[idx], sasa_norm[idx]),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=7, fontweight='bold', color='black',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7, edgecolor='none'))

    ax5.set_xlabel('Normalized Rg', fontsize=11)
    ax5.set_ylabel('Normalized SASA', fontsize=11)
    ax5.set_title('Extreme Conformers (with IDs)', fontsize=12, fontweight='bold')
    ax5.legend(fontsize=8, loc='best')
    ax5.grid(True, alpha=0.3)

    # 6. Time evolution
    ax6 = plt.subplot(2, 3, 6)
    frames = np.arange(len(rg))
    ax6.scatter(frames, rg, c=sasa, s=5, cmap='viridis', alpha=0.5)
    ax6.scatter(selected_indices, rg[selected_indices], c='red', s=30,
               marker='x', label='Selected')
    ax6.set_xlabel('Frame', fontsize=11)
    ax6.set_ylabel('Rg (nm)', fontsize=11)
    ax6.set_title('Rg vs Time (colored by SASA)', fontsize=12, fontweight='bold')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    cbar = plt.colorbar(ax6.collections[0], ax=ax6)
    cbar.set_label('SASA (nm²)', fontsize=10)

    plt.tight_layout()
    plot_file = output_file.replace('.dat', '_plots.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"✓ Plots saved to {plot_file}")
    plt.close()


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Conformational sampling using Rg-SASA grid",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (10×10 grid)
  python conformational_sampling.py npt.xtc npt.gro

  # Custom grid size with plotting
  python conformational_sampling.py npt.xtc npt.gro --grid-size 15 --plot

  # Adjust threshold for sparse sampling
  python conformational_sampling.py npt.xtc npt.gro --threshold 0.05

  # Full customization
  python conformational_sampling.py npt.xtc npt.gro -g 12 -t 0.04 \
      -o my_conformers.dat --plot
        """
    )

    parser.add_argument('trajectory',
                       help='Trajectory file (.xtc, .trr)')
    parser.add_argument('topology',
                       help='Topology file (.gro, .pdb, .tpr)')
    parser.add_argument('-g', '--grid-size',
                       type=int,
                       default=10,
                       help='Grid points per dimension (default: 10)')
    parser.add_argument('-t', '--threshold',
                       type=float,
                       default=0.03,
                       help='Distance threshold in normalized space (default: 0.03)')
    parser.add_argument('-o', '--output',
                       default='conformers.dat',
                       help='Output file (default: conformers.dat)')
    parser.add_argument('--plot',
                       action='store_true',
                       help='Generate plots (requires matplotlib)')

    args = parser.parse_args()

    # Validate inputs
    if not Path(args.trajectory).exists():
        print(f"ERROR: Trajectory not found: {args.trajectory}")
        sys.exit(1)

    if not Path(args.topology).exists():
        print(f"ERROR: Topology not found: {args.topology}")
        sys.exit(1)

    if args.grid_size < 2:
        print(f"ERROR: Grid size must be >= 2, got {args.grid_size}")
        sys.exit(1)

    if args.threshold <= 0:
        print(f"ERROR: Threshold must be positive, got {args.threshold}")
        sys.exit(1)

    # Run analysis
    try:
        selected = analyze_conformers(
            args.trajectory,
            args.topology,
            args.grid_size,
            args.threshold,
            args.output,
            args.plot
        )
        print(f"\n✓ Analysis complete! Selected {len(selected)} conformers\n")
    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
