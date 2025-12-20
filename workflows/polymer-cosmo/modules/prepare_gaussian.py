#!/usr/bin/env python3
"""
Prepare Gaussian input files for COSMO-RS calculations

Takes conformer .gjf files (created manually from PDBs in GaussView) and
injects proper Gaussian commands from template for:
1. Gas phase geometry optimization
2. CPCM optimization in water
3. COSMO-RS single point calculation

Usage:
    python prepare_gaussian.py <conformer_dir> [--template job_example.gjf]

Example:
    python prepare_gaussian.py conformers/
    python prepare_gaussian.py conformers/ --template my_template.gjf
"""

import os
import sys
import argparse
from pathlib import Path
import shutil


def prepare_gaussian_inputs(conformer_dir, template_file, backup=True):
    """
    Prepare Gaussian input files from template.

    Parameters:
    -----------
    conformer_dir : str
        Directory containing conformer .gjf files (from GaussView)
    template_file : str
        Path to job_example.gjf template file
    backup : bool
        Create backup of original .gjf files

    Returns:
    --------
    processed_files : list
        List of processed .gjf files
    """

    conformer_dir = Path(conformer_dir)
    template_file = Path(template_file)

    if not conformer_dir.exists():
        print(f"ERROR: Directory not found: {conformer_dir}")
        sys.exit(1)

    if not template_file.exists():
        print(f"ERROR: Template not found: {template_file}")
        sys.exit(1)

    print("=" * 70)
    print("GAUSSIAN INPUT PREPARATION FOR COSMO-RS")
    print("=" * 70)
    print(f"Conformer directory: {conformer_dir}")
    print(f"Template file:       {template_file}")
    print("=" * 70)

    # Read template
    print(f"\nReading template: {template_file}")
    try:
        with open(template_file, 'r') as f:
            template_content = f.readlines()
    except Exception as e:
        print(f"ERROR: Failed to read template: {e}")
        sys.exit(1)

    print(f"✓ Template loaded ({len(template_content)} lines)")

    # Find all .gjf files in directory
    gjf_files = sorted(conformer_dir.glob("*.gjf"))

    if len(gjf_files) == 0:
        print(f"\nERROR: No .gjf files found in {conformer_dir}")
        print("  Make sure you've converted PDB files to GJF using GaussView first!")
        sys.exit(1)

    print(f"\nFound {len(gjf_files)} .gjf files")

    # Create backup directory if requested
    if backup:
        backup_dir = conformer_dir / "backup_original_gjf"
        backup_dir.mkdir(exist_ok=True)
        print(f"\nCreating backups in: {backup_dir}")

    # Process each .gjf file
    processed_files = []
    print(f"\nProcessing files...")
    print("-" * 70)

    for gjf_file in gjf_files:
        try:
            # Backup original
            if backup:
                shutil.copy2(gjf_file, backup_dir / gjf_file.name)

            # Read conformer structure (skip first 5 lines, take rest)
            with open(gjf_file, 'r') as f:
                # Skip first 5 lines (Gaussian header from GaussView)
                for _ in range(5):
                    f.readline()
                # Read molecular structure
                molecule_structure = f.read()

            # Create new content from template
            new_content = template_content.copy()

            # Replace placeholder name (A7) with conformer name
            conformer_name = gjf_file.stem  # e.g., "conformer_001_frame23"

            for i in range(len(new_content)):
                new_content[i] = new_content[i].replace('A7', conformer_name)

            # Replace *structure* placeholder with actual structure
            # Find line with *structure* and replace
            for i in range(len(new_content)):
                if '*structure*' in new_content[i]:
                    new_content[i] = molecule_structure
                    break

            # Write updated file
            with open(gjf_file, 'w') as f:
                f.writelines(new_content)

            print(f"  ✓ {gjf_file.name:40s} → {conformer_name}")
            processed_files.append(gjf_file)

        except Exception as e:
            print(f"  ✗ {gjf_file.name:40s} ERROR: {e}")
            continue

    print("-" * 70)
    print(f"\n✓ Processed {len(processed_files)} / {len(gjf_files)} files")

    if backup:
        print(f"✓ Original files backed up to: {backup_dir}")

    print("\n" + "=" * 70)
    print("GAUSSIAN FILES READY FOR SUBMISSION")
    print("=" * 70)
    print("\nEach .gjf file now contains:")
    print("  1. Gas phase geometry optimization (BP86/TZVP)")
    print("  2. CPCM optimization in water")
    print("  3. COSMO-RS single point calculation")
    print("\nNext step - generate submission script:")
    print("  python gaussian_submission_generator.py conformers/")
    print("\nOr use your cluster's job scheduler manually.")
    print("=" * 70)

    return processed_files


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Prepare Gaussian input files for COSMO-RS",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Workflow:
---------
1. Extract conformers from MD trajectory using conformational_sampling.py
2. Convert PDB files to .gjf using GaussView (manually)
3. Run this script to add COSMO-RS commands
4. Submit .gjf files to Gaussian

Example:
--------
# After extracting conformers and creating .gjf files in GaussView:
python prepare_gaussian.py conformers/

# Use custom template:
python prepare_gaussian.py conformers/ --template my_template.gjf

# Disable backup:
python prepare_gaussian.py conformers/ --no-backup

Important:
----------
- Original .gjf files should be created in GaussView from PDB files
- This script MODIFIES files in-place (backups created by default)
- Template must have *structure* placeholder
- Template name references (A7) will be replaced with conformer names
        """
    )

    parser.add_argument('conformer_dir',
                       help='Directory containing conformer .gjf files')
    parser.add_argument('-t', '--template',
                       default='../templates/job_example.gjf',
                       help='Path to template file (default: ../templates/job_example.gjf)')
    parser.add_argument('--no-backup',
                       action='store_true',
                       help='Do not create backup of original files')

    args = parser.parse_args()

    # Validate inputs
    conformer_dir = Path(args.conformer_dir)
    if not conformer_dir.exists():
        print(f"ERROR: Directory not found: {conformer_dir}")
        sys.exit(1)

    # Find template (try relative to script location if not found)
    template_file = Path(args.template)
    if not template_file.exists():
        # Try relative to this script
        script_dir = Path(__file__).parent
        template_file = script_dir.parent / "templates" / "job_example.gjf"

    if not template_file.exists():
        print(f"ERROR: Template not found: {args.template}")
        print(f"  Also tried: {template_file}")
        print("\nMake sure job_example.gjf exists in workflows/polymer-cosmo/templates/")
        sys.exit(1)

    # Process files
    try:
        processed = prepare_gaussian_inputs(
            conformer_dir,
            template_file,
            backup=not args.no_backup
        )

        if len(processed) == 0:
            print("\n⚠️  No files were processed!")
            sys.exit(1)

        print(f"\n✓ Successfully prepared {len(processed)} Gaussian input files\n")

    except KeyboardInterrupt:
        print("\n\nInterrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
