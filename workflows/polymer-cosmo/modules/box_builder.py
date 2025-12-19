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

    def _get_atoms_per_molecule(self, gro_file: Path) -> int:
        """Get number of atoms in a single molecule from .gro file."""
        with open(gro_file) as f:
            lines = f.readlines()
            if len(lines) < 2:
                return 0
            try:
                return int(lines[1].strip())
            except ValueError:
                return 0

    def _count_atoms_in_gro(self, gro_file: Path) -> int:
        """Count total atoms in .gro file."""
        with open(gro_file) as f:
            lines = f.readlines()
            if len(lines) < 2:
                return 0
            try:
                return int(lines[1].strip())
            except ValueError:
                return 0

    def add_solvent(self, boxed_gro: Path, n_molecules: int = 15000,
                    timeout_minutes: float = 2.0) -> Tuple[Path, int]:
        """
        Add solvent molecules to box with intelligent retry logic.

        Args:
            boxed_gro: Path to boxed polymer structure
            n_molecules: Target number of solvent molecules to add
            timeout_minutes: Timeout for insertion attempt (default: 2 min)

        Returns:
            Tuple of (Path to solvated structure, actual number of molecules added)
        """
        import time
        import signal
        from contextlib import contextmanager

        solvent_gro = self.work_dir / "solvent.acpype" / "solvent_GMX.gro"
        if not solvent_gro.exists():
            raise FileNotFoundError(f"Solvent structure not found: {solvent_gro}")

        solvated_gro = self.work_dir / "system_solvated.gro"

        # Get atoms per solvent molecule
        atoms_per_solvent = self._get_atoms_per_molecule(solvent_gro)
        initial_atoms = self._count_atoms_in_gro(boxed_gro)

        print(f"\nSolvent molecule has {atoms_per_solvent} atoms")
        print(f"Initial system has {initial_atoms} atoms")

        # Timeout handler
        @contextmanager
        def timeout_context(seconds):
            def timeout_handler(signum, frame):
                raise TimeoutError()

            old_handler = signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(int(seconds))
            try:
                yield
            finally:
                signal.alarm(0)
                signal.signal(signal.SIGALRM, old_handler)

        # Iterative insertion with timeout and reduction
        target = n_molecules
        max_attempts = 5
        timeout_sec = int(timeout_minutes * 60)

        for attempt in range(1, max_attempts + 1):
            print(f"\n--- Attempt {attempt}: Trying to insert {target} molecules ---")
            print(f"Timeout: {timeout_minutes} minutes, Max tries per molecule: 100")

            cmd = (
                f"gmx insert-molecules -f {boxed_gro} -ci {solvent_gro} "
                f"-nmol {target} -o {solvated_gro} -try 100"
            )

            try:
                with timeout_context(timeout_sec):
                    returncode, stdout, stderr = self.run_command(
                        cmd,
                        f"Adding {target} solvent molecules (attempt {attempt})"
                    )

                    if returncode == 0:
                        # Success! Count molecules
                        final_atoms = self._count_atoms_in_gro(solvated_gro)
                        added_atoms = final_atoms - initial_atoms
                        actual_molecules = added_atoms // atoms_per_solvent

                        print(f"\n✓ Successfully inserted solvent molecules")
                        print(f"  Target: {target}")
                        print(f"  Actual molecules added: {actual_molecules}")
                        print(f"  Total atoms: {final_atoms}")
                        print(f"  Fill efficiency: {actual_molecules/target*100:.1f}%")

                        return solvated_gro, actual_molecules

            except TimeoutError:
                print(f"\n⏱ Timeout after {timeout_minutes} minutes")

                # Check if any molecules were added
                if solvated_gro.exists():
                    final_atoms = self._count_atoms_in_gro(solvated_gro)
                    added_atoms = final_atoms - initial_atoms
                    actual_molecules = added_atoms // atoms_per_solvent

                    if actual_molecules > 0:
                        print(f"  Partial success: {actual_molecules} molecules added so far")

                        # If we got >80% of target, accept it
                        if actual_molecules / target > 0.8:
                            print(f"  ✓ Accepting {actual_molecules} molecules (>80% of target)")
                            return solvated_gro, actual_molecules

            except Exception as e:
                print(f"\n✗ Error during insertion: {e}")

            # Reduce target for next attempt
            target = int(target * 0.7)
            print(f"  → Reducing target to {target} molecules for next attempt")

            if target < 100:
                print(f"\n✗ Target too low ({target} < 100). Stopping.")
                break

        # If we get here, all attempts failed
        raise RuntimeError(f"Failed to insert solvent after {max_attempts} attempts")

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

    def _merge_atomtypes(self, polymer_atomtypes: Path, solvent_atomtypes: Path) -> Path:
        """
        Merge atomtypes from polymer and solvent, removing duplicates.

        Args:
            polymer_atomtypes: Path to polymer atomtypes.itp
            solvent_atomtypes: Path to solvent atomtypes.itp

        Returns:
            Path to merged atomtypes file
        """
        merged_file = self.work_dir / "atomtypes.itp"

        # Read atomtypes from both files
        atomtypes = {}  # Dictionary to track unique atomtypes

        def parse_atomtypes(file_path):
            """Parse atomtypes from file."""
            types = {}
            in_atomtypes = False

            with open(file_path) as f:
                for line in f:
                    stripped = line.strip()

                    if stripped.startswith('[') and 'atomtypes' in stripped:
                        in_atomtypes = True
                        continue

                    if in_atomtypes:
                        # Stop at next section
                        if stripped.startswith('['):
                            break

                        # Skip comments and empty lines
                        if not stripped or stripped.startswith(';'):
                            continue

                        # Parse atomtype line
                        parts = line.split()
                        if len(parts) >= 2:
                            atomtype_name = parts[0]
                            types[atomtype_name] = line.rstrip()

            return types

        # Parse both files
        atomtypes.update(parse_atomtypes(polymer_atomtypes))
        atomtypes.update(parse_atomtypes(solvent_atomtypes))

        # Write merged atomtypes with [ defaults ] section
        with open(merged_file, 'w') as f:
            # Add defaults section (required for forcefield)
            f.write("; Combined atomtypes for polymer-solvent system\n")
            f.write("; Generated by polymer-cosmo workflow\n\n")
            f.write("[ defaults ]\n")
            f.write("; nbfunc\tcomb-rule\tgen-pairs\tfudgeLJ\tfudgeQQ\n")
            f.write("1\t2\tyes\t1.000000\t1.000000\n\n")

            # Add atomtypes section
            f.write("[ atomtypes ]\n")
            f.write(";name   bond_type     mass     charge   ptype   sigma         epsilon       Amb\n")

            # Write unique atomtypes sorted by name
            for atomtype_name in sorted(atomtypes.keys()):
                f.write(atomtypes[atomtype_name] + '\n')

        print(f"\n✓ Merged {len(atomtypes)} unique atomtypes into {merged_file.name}")
        return merged_file

    def create_topology(self, polymer_resname: str, solvent_resname: str,
                        n_solvent: int) -> Path:
        """
        Create combined topology file with merged atomtypes.

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

        # Merge atomtypes
        merged_atomtypes = self._merge_atomtypes(polymer_atomtypes, solvent_atomtypes)

        # Check if we need to rename residues in itp files
        polymer_resname_actual = self._get_resname_from_itp(polymer_itp)
        solvent_resname_actual = self._get_resname_from_itp(solvent_itp)

        topology_content = f"""; Topology for polymer-solvent system
; Generated by polymer-cosmo workflow

; Merged atom types (includes [ defaults ])
#include "{merged_atomtypes.name}"

; Polymer topology
#include "polymer.acpype/polymer_GMX.itp"

; Solvent topology
#include "solvent.acpype/solvent_GMX.itp"

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
        print(f"  Polymer: {polymer_resname_actual} (1 molecule)")
        print(f"  Solvent: {solvent_resname_actual} ({n_solvent} molecules)")
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
        solvated_gro, actual_n_solvent = builder.add_solvent(boxed_gro, args.n_solvent)

        # Create topology
        print("\n" + "=" * 60)
        print("STEP 5: Creating topology")
        print("=" * 60)
        topol = builder.create_topology(
            args.polymer_resname,
            args.solvent_resname,
            actual_n_solvent  # Use actual count from insertion
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
