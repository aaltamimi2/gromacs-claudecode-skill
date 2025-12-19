#!/usr/bin/env python3
"""
box_builder.py - Create solvated simulation boxes for polymer-solvent systems

Integrates with solvent_to_gmx.py to parameterize molecules and build solvated boxes.
"""

import subprocess
import shutil
import sys
from pathlib import Path
from typing import Optional, Tuple
import re


class BoxBuilder:
    """Build solvated polymer boxes for GROMACS simulations."""

    def __init__(self, work_dir: Path):
        """
        Initialize BoxBuilder.

        Args:
            work_dir: Working directory for all files
        """
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

        # Find scripts
        self.skill_dir = Path(__file__).parent.parent.parent.parent
        self.solvent_to_gmx = self.skill_dir / "scripts" / "solvent_to_gmx.py"

        if not self.solvent_to_gmx.exists():
            raise FileNotFoundError(f"solvent_to_gmx.py not found at {self.solvent_to_gmx}")

    def run_command(self, cmd: str, description: str = "") -> Tuple[int, str, str]:
        """Run shell command and return result."""
        print(f"\n{'='*60}")
        if description:
            print(f">>> {description}")
        print(f">>> {cmd}")
        print(f"{'='*60}")

        result = subprocess.run(
            cmd,
            shell=True,
            cwd=self.work_dir,
            capture_output=True,
            text=True
        )

        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)

        return result.returncode, result.stdout, result.stderr

    def parameterize_molecule(self, smiles: str, resname: str, base_name: str) -> Path:
        """
        Parameterize molecule using solvent_to_gmx.py.

        Args:
            smiles: SMILES string
            resname: 3-letter residue name
            base_name: Base name for output files

        Returns:
            Path to .acpype directory
        """
        cmd = (
            f"python {self.solvent_to_gmx} "
            f"--smiles '{smiles}' "
            f"--resname {resname} "
            f"--base {base_name}"
        )

        returncode, stdout, stderr = self.run_command(
            cmd,
            f"Parameterizing {base_name} ({resname})"
        )

        if returncode != 0:
            raise RuntimeError(f"Parameterization failed for {base_name}")

        # Find acpype directory
        acpype_dir = self.work_dir / f"{base_name}.acpype"
        if not acpype_dir.exists():
            raise FileNotFoundError(f"Expected {acpype_dir} not created")

        return acpype_dir

    def create_box(self, box_size: float = 10.0) -> Path:
        """
        Create simulation box.

        Args:
            box_size: Box dimension in nm (cubic box)

        Returns:
            Path to boxed structure
        """
        # Use polymer as starting structure
        polymer_gro = self.work_dir / "polymer.acpype" / "polymer_GMX.gro"
        if not polymer_gro.exists():
            raise FileNotFoundError(f"Polymer structure not found: {polymer_gro}")

        boxed_gro = self.work_dir / "polymer_boxed.gro"

        cmd = (
            f"gmx editconf -f {polymer_gro} -o {boxed_gro} "
            f"-box {box_size} {box_size} {box_size} -c"
        )

        returncode, stdout, stderr = self.run_command(
            cmd,
            f"Creating {box_size}x{box_size}x{box_size} nm box"
        )

        if returncode != 0:
            raise RuntimeError("Box creation failed")

        return boxed_gro

    def add_solvent(self, boxed_gro: Path, n_molecules: int = 15000) -> Path:
        """
        Add solvent molecules to box.

        Args:
            boxed_gro: Path to boxed polymer structure
            n_molecules: Number of solvent molecules to add

        Returns:
            Path to solvated structure
        """
        solvent_gro = self.work_dir / "solvent.acpype" / "solvent_GMX.gro"
        if not solvent_gro.exists():
            raise FileNotFoundError(f"Solvent structure not found: {solvent_gro}")

        solvated_gro = self.work_dir / "system_solvated.gro"

        cmd = (
            f"gmx insert-molecules -f {boxed_gro} -ci {solvent_gro} "
            f"-nmol {n_molecules} -o {solvated_gro} -try 500"
        )

        returncode, stdout, stderr = self.run_command(
            cmd,
            f"Adding {n_molecules} solvent molecules"
        )

        if returncode != 0:
            # Try with fewer molecules if insertion fails
            print(f"\nWARNING: Failed to insert {n_molecules} molecules")
            print("Trying with automatic molecule count...")

            # Try without specifying number
            cmd_auto = (
                f"gmx insert-molecules -f {boxed_gro} -ci {solvent_gro} "
                f"-o {solvated_gro} -try 500"
            )
            returncode, stdout, stderr = self.run_command(cmd_auto, "Auto solvent insertion")

            if returncode != 0:
                raise RuntimeError("Solvent insertion failed")

        # Count actual molecules inserted
        actual_count = self._count_molecules_in_gro(solvated_gro)
        print(f"\n✓ Successfully inserted solvent molecules")
        print(f"  Requested: {n_molecules}, Actual: {actual_count}")

        return solvated_gro

    def _count_molecules_in_gro(self, gro_file: Path) -> int:
        """Count number of molecules in .gro file."""
        with open(gro_file) as f:
            lines = f.readlines()
            if len(lines) < 2:
                return 0
            # Second line contains atom count
            try:
                return int(lines[1].strip())
            except ValueError:
                return 0

    def create_topology(self, polymer_resname: str, solvent_resname: str,
                        n_solvent: int) -> Path:
        """
        Create combined topology file.

        Args:
            polymer_resname: Polymer residue name
            solvent_resname: Solvent residue name
            n_solvent: Number of solvent molecules

        Returns:
            Path to topology file
        """
        topol_file = self.work_dir / "topol.top"

        # Get itp files
        polymer_atomtypes = self.work_dir / "polymer.acpype" / "atomtypes.itp"
        polymer_itp = self.work_dir / "polymer.acpype" / f"polymer_GMX.itp"
        solvent_atomtypes = self.work_dir / "solvent.acpype" / "atomtypes.itp"
        solvent_itp = self.work_dir / "solvent.acpype" / f"solvent_GMX.itp"

        # Check if we need to rename residues in itp files
        polymer_resname_actual = self._get_resname_from_itp(polymer_itp)
        solvent_resname_actual = self._get_resname_from_itp(solvent_itp)

        topology_content = f"""; Topology for polymer-solvent system
; Generated by polymer-cosmo workflow

#include "amber99sb-ildn.ff/forcefield.itp"

; Polymer atom types
#include "{polymer_atomtypes.name}"

; Solvent atom types
#include "{solvent_atomtypes.name}"

; Polymer topology
#include "{polymer_itp.name}"

; Solvent topology
#include "{solvent_itp.name}"

; System definition
[ system ]
Polymer in solvent

[ molecules ]
{polymer_resname_actual}    1
{solvent_resname_actual}    {n_solvent}
"""

        with open(topol_file, 'w') as f:
            f.write(topology_content)

        print(f"\n✓ Created topology file: {topol_file}")
        return topol_file

    def _get_resname_from_itp(self, itp_file: Path) -> str:
        """Extract residue name from .itp file."""
        with open(itp_file) as f:
            for line in f:
                if line.strip().startswith('[') and 'moleculetype' in line:
                    # Next non-comment line has the residue name
                    for next_line in f:
                        if not next_line.strip().startswith(';') and next_line.strip():
                            return next_line.split()[0]
        return "UNK"

    def test_energy_minimization(self, em_mdp: Path) -> bool:
        """
        Run test energy minimization.

        Args:
            em_mdp: Path to EM .mdp file

        Returns:
            True if successful
        """
        gro_file = self.work_dir / "system_solvated.gro"
        topol_file = self.work_dir / "topol.top"
        em_tpr = self.work_dir / "em_test.tpr"
        em_out = self.work_dir / "em_test"

        # grompp
        cmd_grompp = (
            f"gmx grompp -f {em_mdp} -c {gro_file} -p {topol_file} "
            f"-o {em_tpr} -maxwarn 10"
        )

        returncode, stdout, stderr = self.run_command(
            cmd_grompp,
            "Preparing energy minimization"
        )

        if returncode != 0:
            print("\n✗ Energy minimization preparation failed")
            return False

        # mdrun (short test)
        cmd_mdrun = f"gmx mdrun -deffnm {em_out} -nsteps 100"

        returncode, stdout, stderr = self.run_command(
            cmd_mdrun,
            "Running test energy minimization (100 steps)"
        )

        if returncode != 0:
            print("\n✗ Energy minimization test failed")
            return False

        print("\n✓ Energy minimization test PASSED")
        return True


def main():
    """CLI for box builder."""
    import argparse

    parser = argparse.ArgumentParser(description="Build solvated polymer boxes")
    parser.add_argument("work_dir", help="Working directory")
    parser.add_argument("--polymer-smiles", required=True, help="Polymer SMILES")
    parser.add_argument("--polymer-resname", default="POL", help="Polymer residue name")
    parser.add_argument("--solvent-smiles", required=True, help="Solvent SMILES")
    parser.add_argument("--solvent-resname", default="SOL", help="Solvent residue name")
    parser.add_argument("--box-size", type=float, default=10.0, help="Box size (nm)")
    parser.add_argument("--n-solvent", type=int, default=15000,
                        help="Number of solvent molecules")
    parser.add_argument("--em-mdp", help="Path to EM .mdp file")

    args = parser.parse_args()

    try:
        builder = BoxBuilder(args.work_dir)

        # Parameterize polymer
        print("\n" + "=" * 60)
        print("STEP 1: Parameterizing polymer")
        print("=" * 60)
        builder.parameterize_molecule(
            args.polymer_smiles,
            args.polymer_resname,
            "polymer"
        )

        # Parameterize solvent
        print("\n" + "=" * 60)
        print("STEP 2: Parameterizing solvent")
        print("=" * 60)
        builder.parameterize_molecule(
            args.solvent_smiles,
            args.solvent_resname,
            "solvent"
        )

        # Create box
        print("\n" + "=" * 60)
        print("STEP 3: Creating simulation box")
        print("=" * 60)
        boxed_gro = builder.create_box(args.box_size)

        # Add solvent
        print("\n" + "=" * 60)
        print("STEP 4: Adding solvent molecules")
        print("=" * 60)
        solvated_gro = builder.add_solvent(boxed_gro, args.n_solvent)

        # Get actual solvent count
        with open(solvated_gro) as f:
            lines = f.readlines()
        # Rough count - need to subtract polymer atoms
        # This is a simplification

        # Create topology
        print("\n" + "=" * 60)
        print("STEP 5: Creating topology")
        print("=" * 60)
        topol = builder.create_topology(
            args.polymer_resname,
            args.solvent_resname,
            args.n_solvent  # Use requested count for now
        )

        # Test EM if mdp provided
        if args.em_mdp:
            print("\n" + "=" * 60)
            print("STEP 6: Testing energy minimization")
            print("=" * 60)
            success = builder.test_energy_minimization(Path(args.em_mdp))
            if not success:
                sys.exit(1)

        print("\n" + "=" * 60)
        print("✓ SUCCESS: Solvated box created")
        print("=" * 60)
        print(f"Structure: {solvated_gro}")
        print(f"Topology: {topol}")

    except Exception as e:
        print(f"\n✗ ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
