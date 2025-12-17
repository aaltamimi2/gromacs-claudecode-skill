#!/usr/bin/env python3
"""
solvent_to_gmx.py - Convert SMILES to GROMACS topology and coordinates

Given a SMILES string, this script:
1. Generates 3D coordinates with RDKit
2. Converts to MOL2 format
3. Runs ACPYPE to generate GROMACS topology (GAFF)
4. Post-processes outputs (normalizes residue names, extracts atomtypes)

Usage:
    python solvent_to_gmx.py --smiles "CCO" --resname ETH --base ethanol
    python solvent_to_gmx.py --smiles "O=C([O-])C1=C(O)C(/N=N/C2=CC=C(C)C=C2S(=O)([O-])=O)=C3C=CC=CC3=C1" --resname DYE

Requirements:
    pip install rdkit openbabel-wheel
    # ACPYPE must be installed and in PATH
"""

import argparse
import subprocess
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


def smiles_to_3d_sdf(smiles: str, out_sdf: Path, title: str, max_attempts: int = 5):
    """
    Generate 3D coordinates from SMILES using RDKit and save as SDF.

    Args:
        smiles: SMILES string
        out_sdf: Output SDF file path
        title: Molecule title for SDF
        max_attempts: Maximum embedding attempts with different random seeds
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"RDKit could not parse SMILES: {smiles}")

    mol = Chem.AddHs(mol)

    # Try embedding with different random seeds if it fails
    success = False
    for attempt in range(max_attempts):
        params = AllChem.ETKDGv3()
        params.randomSeed = 12345 + attempt

        if AllChem.EmbedMolecule(mol, params) == 0:
            success = True
            break

        if attempt < max_attempts - 1:
            print(f"[RDKit] Embedding attempt {attempt + 1} failed, retrying...")

    if not success:
        raise RuntimeError(
            f"RDKit embedding failed after {max_attempts} attempts. "
            "Try simplifying the SMILES or checking for invalid stereochemistry."
        )

    # Optimize geometry
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    else:
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)

    mol.SetProp("_Name", title)

    w = Chem.SDWriter(str(out_sdf))
    w.write(mol)
    w.close()


def run(cmd: list[str]):
    """Execute shell command and raise on error."""
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if p.returncode != 0:
        sys.stderr.write(p.stdout)
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return p.stdout


def normalize_acpype_resname(resname: str):
    """
    Replace UNL -> resname in <resname>.acpype/<resname>_GMX.itp and <resname>_GMX.gro
    """
    acpype_dir = Path(f"{resname}.acpype")
    if not acpype_dir.is_dir():
        raise RuntimeError(f"ACPYPE output directory not found: {acpype_dir}")

    itp = acpype_dir / f"{resname}_GMX.itp"
    gro = acpype_dir / f"{resname}_GMX.gro"

    for f in (itp, gro):
        if not f.exists():
            raise RuntimeError(f"Expected ACPYPE file not found: {f}")

        text = f.read_text()
        text = text.replace("UNL", resname)
        f.write_text(text)

    print(f"[POST] Normalized residue name UNL → {resname} in ACPYPE outputs: {acpype_dir}")


def _find_section_bounds(lines: list[str], section: str) -> tuple[int, int] | None:
    """
    Return (start_idx, end_idx) for a section like 'atomtypes' in a GROMACS .itp.
    """
    header = f"[ {section} ]"
    start = None
    for i, ln in enumerate(lines):
        if ln.strip() == header:
            start = i
            break
    if start is None:
        return None

    end = len(lines)
    for j in range(start + 1, len(lines)):
        s = lines[j].strip()
        if s.startswith("[") and s.endswith("]"):
            end = j
            break
    return start, end


def extract_atomtypes_to_file(resname: str):
    """
    Extract [ atomtypes ] section to separate file and clean up main .itp
    """
    acpype_dir = Path(f"{resname}.acpype")
    itp = acpype_dir / f"{resname}_GMX.itp"
    if not itp.exists():
        raise RuntimeError(f"Expected ACPYPE file not found: {itp}")

    backup = acpype_dir / f"{resname}-original.itp"
    atomtypes_out = acpype_dir / "atomtypes.itp"

    original_text = itp.read_text()
    backup.write_text(original_text)

    lines = original_text.splitlines(keepends=True)

    bounds = _find_section_bounds(lines, "atomtypes")
    if bounds is None:
        raise RuntimeError(f"No [ atomtypes ] section found in {itp}")

    start, end = bounds

    atomtypes_lines = lines[start:end]
    remaining_lines = lines[:start] + lines[end:]

    # Write atomtypes.itp
    atomtypes_out.write_text("".join(atomtypes_lines).rstrip() + "\n")

    # Write cleaned *_GMX.itp (no atomtypes block)
    itp.write_text("".join(remaining_lines).rstrip() + "\n")

    print(f"[POST] Wrote atomtypes to: {atomtypes_out}")
    print(f"[POST] Backed up original itp to: {backup}")
    print(f"[POST] Removed [ atomtypes ] from: {itp}")


def calculate_net_charge(smiles: str) -> int:
    """Calculate net charge from SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    return Chem.GetFormalCharge(mol)


def main():
    ap = argparse.ArgumentParser(
        description="Convert SMILES to GROMACS topology and coordinates using RDKit + ACPYPE"
    )
    ap.add_argument("--smiles", required=True, help="SMILES string of the molecule")
    ap.add_argument("--base", default="MOL", help="Base name for outputs (default: MOL)")
    ap.add_argument("--resname", default="MOL", help="Residue name for GROMACS (default: MOL)")
    ap.add_argument("--title", help="Molecule title (default: same as resname)")
    args = ap.parse_args()

    title = args.title if args.title else args.resname
    base = Path(args.base)
    out_sdf = base.with_suffix(".sdf")
    out_mol2 = base.with_suffix(".mol2")

    print(f"=== SMILES to GROMACS: {title} ===")
    print(f"SMILES: {args.smiles}")

    # Calculate net charge
    net_charge = calculate_net_charge(args.smiles)
    print(f"[RDKit] Detected net charge: {net_charge:+d}")

    # Generate 3D structure
    smiles_to_3d_sdf(args.smiles, out_sdf, title=title)
    print(f"[RDKit] Wrote 3D SDF: {out_sdf}")

    # Convert to MOL2
    run(["obabel", str(out_sdf), "-O", str(out_mol2)])
    print(f"[OpenBabel] Wrote MOL2: {out_mol2}")

    # Run ACPYPE
    acpype_cmd = ["acpype", "-i", str(out_mol2), "-b", args.resname, "-o", "gmx"]
    if net_charge != 0:
        acpype_cmd.extend(["-n", str(net_charge)])
        print(f"[ACPYPE] Using net charge: {net_charge:+d}")

    run(acpype_cmd)
    print(f"[ACPYPE] Done. See: {args.resname}.acpype/")

    # Post-process
    normalize_acpype_resname(args.resname)
    extract_atomtypes_to_file(args.resname)

    print(f"\n✓ SMILES to GROMACS conversion complete!")
    print(f"  SMILES: {args.smiles}")
    print(f"  Net charge: {net_charge:+d}")
    print(f"  Topology: {args.resname}.acpype/{args.resname}_GMX.itp")
    print(f"  Coordinates: {args.resname}.acpype/{args.resname}_GMX.gro")
    print(f"  Atomtypes: {args.resname}.acpype/atomtypes.itp")
    print(f"\nNext steps:")
    print(f"  1. Include atomtypes in your topology: #include \"{args.resname}.acpype/atomtypes.itp\"")
    print(f"  2. Include molecule: #include \"{args.resname}.acpype/{args.resname}_GMX.itp\"")
    print(f"  3. Add to [ molecules ]: {args.resname}  N")

    if net_charge != 0:
        print(f"\n⚠ WARNING: Molecule has net charge {net_charge:+d}")
        print(f"  You may need to add counterions for charge neutralization")


if __name__ == "__main__":
    main()
