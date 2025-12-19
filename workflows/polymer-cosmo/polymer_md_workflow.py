#!/usr/bin/env python3
"""
polymer_md_workflow.py - Main orchestrator for polymer-solvent MD simulations

End-to-end workflow:
1. Generate polymer oligomer (validate backbone atoms)
2. Validate/get solvent SMILES
3. Parameterize polymer and solvent
4. Create solvated simulation box
5. Test energy minimization
6. Generate submission scripts (swarm/workstation)

Usage:
    python polymer_md_workflow.py --polymer PVDF --solvent NMP --work-dir ./pvdf_nmp
"""

import sys
import shutil
from pathlib import Path
from typing import Optional, Tuple
import argparse

# Add modules to path
WORKFLOW_DIR = Path(__file__).parent
sys.path.insert(0, str(WORKFLOW_DIR / "modules"))

from polymer_builder import (
    create_polymer_structure,
    recommend_solvents,
    validate_solvent_smiles,
    get_polymer_smiles,
    POLYMER_DATABASE,
    SOLVENT_DATABASE
)
from box_builder import BoxBuilder
from submission_generator import SubmissionGenerator


class PolymerMDWorkflow:
    """Main workflow orchestrator for polymer MD simulations."""

    def __init__(self, work_dir: Path, interactive: bool = True):
        """
        Initialize workflow.

        Args:
            work_dir: Working directory for all files
            interactive: If True, prompt for user confirmation at checkpoints
        """
        self.work_dir = Path(work_dir)
        self.interactive = interactive
        self.template_dir = WORKFLOW_DIR / "templates"

        # State tracking
        self.polymer_info = None
        self.solvent_info = None

    def _print_section(self, title: str):
        """Print formatted section header."""
        print("\n" + "=" * 70)
        print(f"  {title}")
        print("=" * 70)

    def _confirm(self, prompt: str = "Continue?") -> bool:
        """Ask user for confirmation."""
        if not self.interactive:
            return True

        response = input(f"\n{prompt} [Y/n]: ").strip().lower()
        return response in ['', 'y', 'yes']

    def step1_generate_polymer(
        self,
        polymer_input: str,
        oligomer_size: Optional[int] = None,
        target_backbone: int = 11,
        max_backbone: int = 20
    ) -> dict:
        """
        Step 1: Generate and validate polymer oligomer.

        Args:
            polymer_input: Polymer name or SMILES
            oligomer_size: Explicit oligomer size (overrides auto)
            target_backbone: Target backbone atoms
            max_backbone: Maximum backbone atoms

        Returns:
            dict with polymer info
        """
        self._print_section("STEP 1: Generate Polymer Oligomer")

        # Create polymer structure
        polymer_info = create_polymer_structure(
            polymer_input,
            oligomer_size,
            target_backbone,
            max_backbone
        )

        print(f"\n✓ Polymer: {polymer_info['polymer_name']}")
        print(f"  Monomer SMILES: {polymer_info['monomer_smiles']}")
        print(f"  Oligomer size: {polymer_info['n_units']} units")
        print(f"  Backbone atoms: {polymer_info['backbone_atoms']} "
              f"(target: {target_backbone}, max: {max_backbone})")
        print(f"  Oligomer SMILES: {polymer_info['oligomer_smiles']}")

        # Recommend solvents
        solvents = recommend_solvents(polymer_info['polymer_name'])
        if solvents:
            print(f"\n  Recommended solvents: {', '.join(solvents)}")

        self.polymer_info = polymer_info
        return polymer_info

    def step2_validate_solvent(self, solvent_input: str) -> Tuple[str, str]:
        """
        Step 2: Validate solvent SMILES.

        Args:
            solvent_input: Solvent name or SMILES

        Returns:
            (solvent_smiles, solvent_name)
        """
        self._print_section("STEP 2: Validate Solvent")

        smiles, message = validate_solvent_smiles(solvent_input)

        if smiles is None:
            print(f"\n✗ {message}")
            raise ValueError(f"Invalid solvent: {solvent_input}")

        print(f"\n✓ {message}")
        print(f"  SMILES: {smiles}")

        # Get solvent name
        solvent_name = solvent_input if solvent_input.lower() in SOLVENT_DATABASE else "CustomSolvent"

        # USER VERIFICATION REQUIRED
        if self.interactive:
            print(f"\n⚠️  IMPORTANT: Verify this SMILES is correct for your solvent!")
            if not self._confirm("Is this solvent SMILES correct?"):
                raise ValueError("Solvent SMILES not confirmed by user")

        self.solvent_info = {"smiles": smiles, "name": solvent_name}
        return smiles, solvent_name

    def step3_build_solvated_box(
        self,
        polymer_smiles: str,
        solvent_smiles: str,
        polymer_resname: str = "POL",
        solvent_resname: str = "SOL",
        box_size: float = 10.0,
        n_solvent: int = 15000
    ):
        """
        Step 3: Build solvated simulation box.

        Args:
            polymer_smiles: Polymer SMILES
            solvent_smiles: Solvent SMILES
            polymer_resname: 3-letter polymer residue name
            solvent_resname: 3-letter solvent residue name
            box_size: Box dimension (nm)
            n_solvent: Number of solvent molecules
        """
        self._print_section("STEP 3: Build Solvated Box")

        builder = BoxBuilder(self.work_dir)

        # Parameterize polymer
        print("\n--- Parameterizing Polymer ---")
        polymer_acpype = builder.parameterize_molecule(
            polymer_smiles,
            polymer_resname,
            "polymer"
        )
        print(f"✓ Polymer parameterization complete: {polymer_acpype}")

        # Parameterize solvent
        print("\n--- Parameterizing Solvent ---")
        solvent_acpype = builder.parameterize_molecule(
            solvent_smiles,
            solvent_resname,
            "solvent"
        )
        print(f"✓ Solvent parameterization complete: {solvent_acpype}")

        # Create box
        print("\n--- Creating Simulation Box ---")
        boxed_gro = builder.create_box(box_size)
        print(f"✓ Box created: {boxed_gro}")

        # Add solvent
        print("\n--- Adding Solvent Molecules ---")
        solvated_gro = builder.add_solvent(boxed_gro, n_solvent)
        print(f"✓ Solvated system: {solvated_gro}")

        # Create topology
        print("\n--- Creating Topology ---")
        topol = builder.create_topology(polymer_resname, solvent_resname, n_solvent)
        print(f"✓ Topology: {topol}")

        if self.interactive:
            if not self._confirm("Box creation complete. Continue?"):
                raise RuntimeError("User stopped workflow")

    def step4_test_energy_minimization(self):
        """Step 4: Test energy minimization."""
        self._print_section("STEP 4: Test Energy Minimization")

        # Copy EM template
        em_template = self.template_dir / "em.mdp"
        em_mdp = self.work_dir / "em.mdp"
        shutil.copy(em_template, em_mdp)
        print(f"✓ Copied EM parameters: {em_mdp}")

        # Run test
        builder = BoxBuilder(self.work_dir)
        success = builder.test_energy_minimization(em_mdp)

        if not success:
            print("\n✗ Energy minimization test FAILED")
            print("  Check topology and structure files")
            if self.interactive:
                if not self._confirm("EM test failed. Continue anyway?"):
                    raise RuntimeError("Energy minimization test failed")
        else:
            print("\n✓ Energy minimization test PASSED")

    def step5_generate_submission_scripts(
        self,
        job_name: str,
        email: str = "aaltamimi2@wisc.edu",
        target: str = "both"
    ):
        """
        Step 5: Generate submission scripts.

        Args:
            job_name: Name for SLURM job
            email: Email for notifications
            target: "swarm", "workstation", or "both"
        """
        self._print_section("STEP 5: Generate Submission Scripts")

        # Copy MDP templates
        for mdp_name in ["nvt.mdp", "npt.mdp"]:
            src = self.template_dir / mdp_name
            dst = self.work_dir / mdp_name
            shutil.copy(src, dst)
            print(f"✓ Copied {mdp_name}")

        # Generate scripts
        gen = SubmissionGenerator(self.work_dir)

        if target in ["swarm", "both"]:
            gen.generate_swarm_script(
                job_name,
                email,
                walltime="72:00:00",
                run_em=True,
                run_nvt=True,
                run_npt=True
            )

        if target in ["workstation", "both"]:
            gen.generate_workstation_script(
                job_name,
                run_em=True,
                run_nvt=True,
                run_npt=True
            )

        gen.generate_preparation_script()

        print("\n✓ All submission scripts generated")

    def run_full_workflow(
        self,
        polymer: str,
        solvent: str,
        oligomer_size: Optional[int] = None,
        job_name: Optional[str] = None,
        email: str = "aaltamimi2@wisc.edu",
        box_size: float = 10.0,
        n_solvent: int = 15000,
        target: str = "both"
    ):
        """
        Run complete workflow from polymer + solvent to submission scripts.

        Args:
            polymer: Polymer name or SMILES
            solvent: Solvent name or SMILES
            oligomer_size: Number of monomer units (auto if None)
            job_name: SLURM job name (auto-generated if None)
            email: Email for notifications
            box_size: Simulation box size (nm)
            n_solvent: Number of solvent molecules
            target: Target platform (swarm/workstation/both)
        """
        try:
            # Step 1: Generate polymer
            polymer_info = self.step1_generate_polymer(polymer, oligomer_size)

            if self.interactive:
                if not self._confirm("Polymer structure OK?"):
                    print("Workflow stopped by user")
                    return

            # Step 2: Validate solvent
            solvent_smiles, solvent_name = self.step2_validate_solvent(solvent)

            # Step 3: Build solvated box
            self.step3_build_solvated_box(
                polymer_info['oligomer_smiles'],
                solvent_smiles,
                polymer_resname="POL",
                solvent_resname="SOL",
                box_size=box_size,
                n_solvent=n_solvent
            )

            # Step 4: Test EM
            self.step4_test_energy_minimization()

            # Step 5: Generate scripts
            if job_name is None:
                job_name = f"{polymer_info['polymer_name']}_{solvent_name}_md"

            self.step5_generate_submission_scripts(job_name, email, target)

            # Success summary
            self._print_section("✓ WORKFLOW COMPLETE")
            print(f"\nWorking directory: {self.work_dir}")
            print(f"\nGenerated files:")
            print(f"  - system_solvated.gro    (solvated structure)")
            print(f"  - topol.top              (topology)")
            print(f"  - *.mdp                  (MD parameters)")
            print(f"  - submit_swarm.sh        (swarm submission)")
            print(f"  - run_workstation.sh     (workstation script)")
            print(f"  - prepare_simulations.sh (prep .tpr files)")

            print(f"\n📋 Next Steps:")
            print(f"  1. Review files in: {self.work_dir}")
            print(f"  2. Prepare .tpr files: ./prepare_simulations.sh")
            print(f"  3. Transfer to swarm: scp * swarm:~/project_dir/")
            print(f"  4. Submit: ssh swarm 'cd ~/project_dir && sbatch submit_swarm.sh'")

        except Exception as e:
            self._print_section("✗ WORKFLOW FAILED")
            print(f"\nError: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Polymer-Solvent MD Workflow",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # PVDF 4mer in NMP (auto oligomer size)
  python polymer_md_workflow.py --polymer PVDF --solvent NMP --work-dir ./pvdf_nmp

  # Custom oligomer size
  python polymer_md_workflow.py --polymer PVDF --solvent NMP --oligomer-size 4 --work-dir ./pvdf4_nmp

  # Custom polymer SMILES
  python polymer_md_workflow.py --polymer "C(C(F)(F))([H])[H]" --solvent "CN1CCCC1=O" --work-dir ./custom

  # Non-interactive mode (no prompts)
  python polymer_md_workflow.py --polymer PVDF --solvent NMP --work-dir ./test --no-interactive

Supported polymers: """ + ", ".join(POLYMER_DATABASE.keys()) + """
Supported solvents: """ + ", ".join(SOLVENT_DATABASE.keys())
    )

    parser.add_argument("--polymer", required=True,
                        help="Polymer name or SMILES")
    parser.add_argument("--solvent", required=True,
                        help="Solvent name or SMILES")
    parser.add_argument("--work-dir", required=True,
                        help="Working directory for simulation files")
    parser.add_argument("--oligomer-size", type=int,
                        help="Number of monomer units (auto if not specified)")
    parser.add_argument("--job-name",
                        help="SLURM job name (auto-generated if not specified)")
    parser.add_argument("--email", default="aaltamimi2@wisc.edu",
                        help="Email for job notifications")
    parser.add_argument("--box-size", type=float, default=10.0,
                        help="Simulation box size in nm (default: 10.0)")
    parser.add_argument("--n-solvent", type=int, default=15000,
                        help="Number of solvent molecules (default: 15000)")
    parser.add_argument("--target", choices=["swarm", "workstation", "both"],
                        default="both",
                        help="Target platform for scripts (default: both)")
    parser.add_argument("--no-interactive", action="store_true",
                        help="Run without user prompts (non-interactive)")

    args = parser.parse_args()

    # Create workflow
    workflow = PolymerMDWorkflow(
        work_dir=args.work_dir,
        interactive=not args.no_interactive
    )

    # Run workflow
    workflow.run_full_workflow(
        polymer=args.polymer,
        solvent=args.solvent,
        oligomer_size=args.oligomer_size,
        job_name=args.job_name,
        email=args.email,
        box_size=args.box_size,
        n_solvent=args.n_solvent,
        target=args.target
    )


if __name__ == "__main__":
    main()
