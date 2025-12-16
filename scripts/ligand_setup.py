#!/usr/bin/env python3
"""
ligand_setup.py - Automated ligand parameterization for GROMACS

Given a compound name, this script:
1. Fetches SMILES from PubChem
2. Generates 3D coordinates with RDKit
3. Converts to MOL2 format
4. Runs ACPYPE to generate GROMACS topology (GAFF)
5. Post-processes outputs (normalizes residue names, extracts atomtypes)

Usage:
    python ligand_setup.py --name "hexane" --resname HEX
    python ligand_setup.py --name "benzene" --base benzene --resname BNZ

Requirements:
    pip install pubchempy rdkit openbabel-wheel
    # ACPYPE must be installed and in PATH
"""

import argparse
import subprocess
import sys
from pathlib import Path

import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem


def name_to_smiles(name: str) -> str:
    """Fetch canonical SMILES from PubChem by compound name."""
    comps = pcp.get_compounds(name, "name")
    if not comps:
        raise RuntimeError(f"No PubChem result for name: {name!r}")
    smi = comps[0].canonical_smiles
    if not smi:
        raise RuntimeError(f"PubChem hit had no canonical_smiles for {name!r}")
    return smi


def smiles_to_3d_sdf(smiles: str, out_sdf: Path, title: str):
    """Generate 3D coordinates from SMILES using RDKit and save as SDF."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"RDKit could not parse SMILES: {smiles}")

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 12345
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError("RDKit embedding failed (try a different SMILES/protonation).")

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
    The slice is inclusive of the section header line and extends until just before
    the next '[ ... ]' header or EOF.
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
    From <resname>.acpype/<resname>_GMX.itp:
      - Save a backup <resname>-original.itp
      - Extract the entire [ atomtypes ] section (including comments directly under it)
        into <resname>.acpype/atomtypes.itp
      - Remove that section from the original <resname>_GMX.itp (so it becomes a clean molecule .itp)
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


def main():
    ap = argparse.ArgumentParser(
        description="Automated ligand parameterization for GROMACS using PubChem + RDKit + ACPYPE"
    )
    ap.add_argument("--name", required=True, help="Common name, e.g. 'hexane'")
    ap.add_argument("--base", default="LIG", help="Base name for outputs (default: LIG)")
    ap.add_argument("--resname", default="LIG", help="Residue name for ACPYPE (default: LIG)")
    args = ap.parse_args()

    base = Path(args.base)
    out_sdf = base.with_suffix(".sdf")
    out_mol2 = base.with_suffix(".mol2")

    print(f"=== Ligand Setup: {args.name} ===")

    smiles = name_to_smiles(args.name)
    print(f"[PubChem] {args.name} -> SMILES: {smiles}")

    smiles_to_3d_sdf(smiles, out_sdf, title=args.name)
    print(f"[RDKit] Wrote 3D SDF: {out_sdf}")

    run(["obabel", str(out_sdf), "-O", str(out_mol2)])
    print(f"[OpenBabel] Wrote MOL2: {out_mol2}")

    run(["acpype", "-i", str(out_mol2), "-b", args.resname, "-o", "gmx"])
    print(f"[ACPYPE] Done. See: {args.resname}.acpype/")

    normalize_acpype_resname(args.resname)
    extract_atomtypes_to_file(args.resname)

    print(f"\n✓ Ligand setup complete!")
    print(f"  Topology: {args.resname}.acpype/{args.resname}_GMX.itp")
    print(f"  Coordinates: {args.resname}.acpype/{args.resname}_GMX.gro")
    print(f"  Atomtypes: {args.resname}.acpype/atomtypes.itp")
    print(f"\nNext steps:")
    print(f"  1. Include atomtypes in your topology: #include \"{args.resname}.acpype/atomtypes.itp\"")
    print(f"  2. Include molecule: #include \"{args.resname}.acpype/{args.resname}_GMX.itp\"")
    print(f"  3. Add to [ molecules ]: {args.resname}  1")


if __name__ == "__main__":
    main()
