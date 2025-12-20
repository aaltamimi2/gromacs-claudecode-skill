#!/usr/bin/env python3
"""
polymer_builder.py - Generate polymer oligomers with backbone atom validation

Creates oligomers from monomer SMILES with backbone atom constraints.
Integrates with solvent_to_gmx.py for GROMACS parameterization.
"""

import subprocess
import sys
from pathlib import Path
from typing import Tuple, Optional

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
except ImportError:
    print("ERROR: RDKit not found. Install with: conda install -c conda-forge rdkit")
    sys.exit(1)


# Common polymer SMILES patterns
POLYMER_DATABASE = {
    "PVDF": {
        "smiles": "C(C(F)(F))([H])[H]",  # -CH2-CF2- monomer
        "name": "Polyvinylidene fluoride",
        "good_solvents": ["NMP", "DMF", "DMSO", "acetone"]
    },
    "PEO": {
        "smiles": "CCO",  # -CH2-CH2-O- monomer
        "name": "Polyethylene oxide",
        "good_solvents": ["water", "chloroform", "acetonitrile"]
    },
    "PMMA": {
        "smiles": "CC(C)(C(=O)OC)",  # Methyl methacrylate monomer
        "name": "Polymethyl methacrylate",
        "good_solvents": ["chloroform", "THF", "toluene"]
    },
    "PS": {
        "smiles": "C(C(c1ccccc1))([H])[H]",  # Styrene monomer
        "name": "Polystyrene",
        "good_solvents": ["toluene", "chloroform", "THF"]
    }
}

SOLVENT_DATABASE = {
    "NMP": "CN1CCCC1=O",  # N-Methyl-2-pyrrolidone
    "DMF": "CN(C)C=O",  # Dimethylformamide
    "DMSO": "CS(=O)C",  # Dimethyl sulfoxide
    "water": "O",
    "chloroform": "ClC(Cl)Cl",
    "acetone": "CC(=O)C",
    "acetonitrile": "CC#N",
    "THF": "C1CCOC1",  # Tetrahydrofuran
    "toluene": "Cc1ccccc1"
}


def count_backbone_atoms(smiles: str) -> int:
    """Count non-hydrogen backbone atoms in a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0

    # Count heavy atoms (non-H)
    return mol.GetNumHeavyAtoms()


def build_oligomer(monomer_smiles: str, n_units: int) -> Tuple[str, int]:
    """
    Build oligomer from monomer SMILES.

    Returns:
        (oligomer_smiles, backbone_atoms)
    """
    # Simple approach: repeat monomer unit
    # For more complex polymers, might need custom linking logic
    mol = Chem.MolFromSmiles(monomer_smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {monomer_smiles}")

    # Create oligomer by combining monomers
    # This is a simplified approach - for real polymers you may need
    # to handle end groups and linking chemistry
    oligomer_smiles = monomer_smiles
    for i in range(1, n_units):
        oligomer_smiles += "." + monomer_smiles

    # Try to create a single molecule (if possible)
    oligomer_mol = Chem.MolFromSmiles(oligomer_smiles)
    if oligomer_mol is None:
        # If can't parse, return as-is
        backbone_atoms = count_backbone_atoms(monomer_smiles) * n_units
        return oligomer_smiles, backbone_atoms

    backbone_atoms = oligomer_mol.GetNumHeavyAtoms()
    return Chem.MolToSmiles(oligomer_mol), backbone_atoms


def recommend_solvents(polymer_name: str) -> list:
    """Recommend good solvents for a polymer."""
    polymer_name = polymer_name.upper()
    if polymer_name in POLYMER_DATABASE:
        return POLYMER_DATABASE[polymer_name]["good_solvents"]
    return []


def validate_solvent_smiles(solvent_name: str) -> Tuple[Optional[str], str]:
    """
    Get or validate solvent SMILES.

    Returns:
        (smiles, validation_message)
    """
    solvent_name_lower = solvent_name.lower()

    if solvent_name_lower in SOLVENT_DATABASE:
        smiles = SOLVENT_DATABASE[solvent_name_lower]
        return smiles, f"Found {solvent_name}: {smiles}"

    # Try to parse as SMILES
    mol = Chem.MolFromSmiles(solvent_name)
    if mol is not None:
        canonical_smiles = Chem.MolToSmiles(mol)
        return canonical_smiles, f"Custom solvent SMILES validated: {canonical_smiles}"

    return None, f"Unknown solvent '{solvent_name}' and not a valid SMILES string"


def get_polymer_smiles(polymer_input: str) -> Tuple[str, str]:
    """
    Get monomer SMILES from polymer name or SMILES input.

    Returns:
        (monomer_smiles, polymer_name)
    """
    polymer_input_upper = polymer_input.upper()

    # Check if it's a known polymer
    if polymer_input_upper in POLYMER_DATABASE:
        data = POLYMER_DATABASE[polymer_input_upper]
        return data["smiles"], polymer_input_upper

    # Try to parse as SMILES
    mol = Chem.MolFromSmiles(polymer_input)
    if mol is not None:
        canonical_smiles = Chem.MolToSmiles(mol)
        return canonical_smiles, "CustomPolymer"

    raise ValueError(f"Unknown polymer '{polymer_input}' and not a valid SMILES")


def calculate_optimal_oligomer_size(monomer_smiles: str,
                                     target_backbone: int = 11,
                                     max_backbone: int = 20) -> int:
    """
    Calculate optimal oligomer size to hit target backbone atoms.

    Args:
        monomer_smiles: SMILES of monomer unit
        target_backbone: Ideal backbone atom count (default 11)
        max_backbone: Maximum allowed backbone atoms (default 20)

    Returns:
        n_units: Number of monomer units
    """
    atoms_per_monomer = count_backbone_atoms(monomer_smiles)
    if atoms_per_monomer == 0:
        raise ValueError("Invalid monomer SMILES")

    # Calculate n for target
    n_target = max(1, round(target_backbone / atoms_per_monomer))

    # Check if it exceeds max
    if n_target * atoms_per_monomer > max_backbone:
        n_target = max(1, int(max_backbone / atoms_per_monomer))

    return n_target


def create_polymer_structure(polymer_name: str,
                              oligomer_size: Optional[int] = None,
                              target_backbone: int = 11,
                              max_backbone: int = 20) -> dict:
    """
    Main function to create polymer oligomer structure.

    Args:
        polymer_name: Name or SMILES of polymer
        oligomer_size: Explicit number of monomer units (overrides auto-calc)
        target_backbone: Target backbone atoms if auto-calculating
        max_backbone: Maximum backbone atoms allowed

    Returns:
        dict with keys: monomer_smiles, oligomer_smiles, n_units, backbone_atoms, polymer_name
    """
    # Get monomer SMILES
    monomer_smiles, name = get_polymer_smiles(polymer_name)

    # Calculate oligomer size if not specified
    if oligomer_size is None:
        oligomer_size = calculate_optimal_oligomer_size(
            monomer_smiles, target_backbone, max_backbone
        )

    # Build oligomer
    oligomer_smiles, backbone_atoms = build_oligomer(monomer_smiles, oligomer_size)

    # Check backbone constraint
    # If oligomer_size=1, user likely provided a complete oligomer SMILES (not a monomer)
    # In this case, warn but don't error
    if backbone_atoms > max_backbone:
        if oligomer_size == 1:
            # User provided complete oligomer - warn but allow
            print(f"\n⚠️  WARNING: Oligomer has {backbone_atoms} backbone atoms (exceeds recommended max {max_backbone})")
            print(f"   This may result in slower simulations and DFT calculations.")
            print(f"   Consider using a smaller oligomer if possible.")
            print(f"   Proceeding anyway...\n")
        else:
            # Auto-generated from monomer - this is an error
            raise ValueError(
                f"Oligomer has {backbone_atoms} backbone atoms, exceeds max {max_backbone}. "
                f"Reduce oligomer size to {int(max_backbone / (backbone_atoms / oligomer_size))} or fewer units."
            )

    return {
        "monomer_smiles": monomer_smiles,
        "oligomer_smiles": oligomer_smiles,
        "n_units": oligomer_size,
        "backbone_atoms": backbone_atoms,
        "polymer_name": name
    }


def main():
    """CLI for testing polymer builder."""
    import argparse

    parser = argparse.ArgumentParser(description="Build polymer oligomers")
    parser.add_argument("polymer", help="Polymer name or SMILES")
    parser.add_argument("-n", "--n-units", type=int, help="Number of monomer units")
    parser.add_argument("--target-backbone", type=int, default=11,
                        help="Target backbone atoms (default: 11)")
    parser.add_argument("--max-backbone", type=int, default=20,
                        help="Max backbone atoms (default: 20)")
    parser.add_argument("--recommend-solvents", action="store_true",
                        help="Show recommended solvents")

    args = parser.parse_args()

    try:
        result = create_polymer_structure(
            args.polymer,
            args.n_units,
            args.target_backbone,
            args.max_backbone
        )

        print(f"\n=== Polymer Oligomer Generated ===")
        print(f"Polymer: {result['polymer_name']}")
        print(f"Monomer SMILES: {result['monomer_smiles']}")
        print(f"Oligomer size: {result['n_units']} units")
        print(f"Backbone atoms: {result['backbone_atoms']}")
        print(f"Oligomer SMILES: {result['oligomer_smiles']}")

        if args.recommend_solvents:
            solvents = recommend_solvents(result['polymer_name'])
            if solvents:
                print(f"\nRecommended solvents: {', '.join(solvents)}")

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
