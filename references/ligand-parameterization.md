# Ligand Parameterization for GROMACS

Guide for generating GROMACS-compatible topology files for small molecules and ligands.

## Quick Start with ligand_setup.py

The repository includes an automated script for ligand parameterization:

```bash
# Basic usage
python scripts/ligand_setup.py --name "hexane" --resname HEX

# Custom base name
python scripts/ligand_setup.py --name "benzene" --base benzene --resname BNZ

# Output:
# - HEX.acpype/HEX_GMX.itp      # Topology
# - HEX.acpype/HEX_GMX.gro      # Coordinates
# - HEX.acpype/atomtypes.itp    # Atomtypes (separate)
```

**What it does:**
1. Fetches structure from PubChem by name
2. Generates 3D coordinates with RDKit
3. Runs ACPYPE to generate GAFF parameters
4. Post-processes to extract atomtypes and normalize residue names

### Requirements

```bash
# Install Python dependencies
pip install pubchempy rdkit openbabel-wheel

# Install ACPYPE
conda install -c conda-forge acpype
# or
pip install acpype
```

### Integrating into GROMACS Topology

```bash
# In your topol.top:
; Include forcefield
#include "amber99sb-ildn.ff/forcefield.itp"

; Include ligand atomtypes BEFORE protein
#include "HEX.acpype/atomtypes.itp"

; Include protein topology
#include "protein.itp"

; Include ligand topology
#include "HEX.acpype/HEX_GMX.itp"

; Include water/ions
#include "amber99sb-ildn.ff/tip3p.itp"

[ system ]
Protein-Ligand Complex

[ molecules ]
Protein_chain_A    1
HEX                1
SOL                5000
NA                 10
CL                 10
```

## Manual Parameterization Methods

### Method 1: ACPYPE (AMBER GAFF)

**Best for:** Small organic molecules, drugs, most ligands

```bash
# 1. Prepare MOL2 file (with correct bond orders)
# Use Avogadro, ChemDraw, or OpenBabel

# 2. Run ACPYPE
acpype -i ligand.mol2 -b LIG -o gmx

# 3. Outputs in LIG.acpype/:
# - LIG_GMX.itp    # Topology
# - LIG_GMX.gro    # Structure
```

**Post-processing:**
```bash
# Normalize residue name (if needed)
cd LIG.acpype/
sed -i 's/UNL/LIG/g' LIG_GMX.itp
sed -i 's/UNL/LIG/g' LIG_GMX.gro

# Extract atomtypes to separate file
# (Manual: copy [ atomtypes ] section to atomtypes.itp)
```

### Method 2: CGenFF (CHARMM)

**Best for:** Use with CHARMM force fields

```bash
# 1. Convert ligand to MOL2
obabel ligand.pdb -O ligand.mol2

# 2. Upload to CGenFF server
# https://cgenff.umaryland.edu/

# 3. Download .str file

# 4. Convert to GROMACS format
python cgenff_charmm2gmx.py LIG ligand.mol2 ligand.str charmm36.ff
```

**Outputs:**
- `LIG.itp` - Topology
- `LIG.prm` - Parameters
- `LIG_ini.pdb` - Structure

### Method 3: ATB (Automated Topology Builder)

**Best for:** GROMOS force fields

```bash
# 1. Upload structure to ATB
# https://atb.uq.edu.au/

# 2. Select GROMOS force field version
# Download .pdb and .itp files

# 3. No post-processing needed (already GROMACS format)
```

### Method 4: LigParGen (OPLS)

**Best for:** OPLS/AA force fields

```bash
# 1. Upload to LigParGen server
# http://zarbi.chem.yale.edu/ligpargen/

# 2. Download GROMACS files
# - ligand.gro
# - ligand.itp

# 3. Include in topology
```

## Protonation States

**Critical:** Ensure correct protonation for target pH (usually 7.0-7.4)

### Using Open Babel
```bash
# Add hydrogens at pH 7.4
obabel input.sdf -O output.sdf -p 7.4
```

### Using RDKit
```python
from rdkit import Chem
from rdkit.Chem import AllChem

mol = Chem.MolFromSmiles("CC(=O)O")  # Acetic acid
mol = Chem.AddHs(mol)  # Add hydrogens

# Or use protomer generation for pH-dependent states
```

### Using ChemAxon
```bash
# Commercial tool, very accurate
cxcalc pka input.sdf
cxcalc majorms -H 7.4 input.sdf
```

## Validation Checklist

After parameterization, verify:

```bash
# 1. Atom count matches
grep -A 100 "atoms" LIG.itp | grep -v "^;" | wc -l
grep ATOM LIG.gro | wc -l

# 2. Residue name consistent
grep "LIG" LIG.itp
grep "LIG" LIG.gro

# 3. Total charge
# Check [ atoms ] section, sum the charges

# 4. Visual inspection
pymol LIG.gro
# Check: bond connectivity, geometry, no clashes

# 5. Test grompp
gmx_mpi grompp -f em.mdp -c LIG.gro -p topol.top -o test.tpr
```

## Common Issues

### Issue: "Atomtype not found"

**Cause:** Atomtypes not included before molecule topology

**Fix:**
```bash
# In topol.top, atomtypes MUST come before molecule
#include "LIG.acpype/atomtypes.itp"   # First
#include "protein.itp"                 # Then protein
#include "LIG.acpype/LIG_GMX.itp"      # Then ligand
```

### Issue: "Non-integer charge"

**Cause:** Ligand has fractional net charge

**Fix:**
1. Check protonation state (add/remove H+)
2. Verify correct oxidation state
3. Add counterion if needed

### Issue: "Bond length constraint failed"

**Cause:** Poor initial geometry or wrong bond orders

**Fix:**
```bash
# Regenerate 3D coordinates
obabel ligand.sdf -O ligand_3d.sdf --gen3d

# Or use RDKit with MMFF optimization
python scripts/ligand_setup.py --name "your_compound"
```

### Issue: "Parameters missing for dihedral"

**Cause:** Uncommon chemical functionality

**Fix:**
1. Use general atom types (edit .itp)
2. Try different parameterization tool
3. Manual parameter derivation (advanced)

## Protein-Ligand Complex Setup

### Full Workflow

```bash
# 1. Parameterize ligand
python scripts/ligand_setup.py --name "ibuprofen" --resname IBU

# 2. Process protein
gmx_mpi pdb2gmx -f protein.pdb -o protein.gro -p topol.top -ff amber99sb-ildn

# 3. Combine structures
# Method A: Manual editing
cat protein.gro > complex.gro
grep -v "Generated" IBU.acpype/IBU_GMX.gro >> complex.gro

# Method B: Use editconf
gmx_mpi editconf -f protein.gro -o protein_box.gro
# Position ligand manually, then combine

# 4. Edit topology
# In topol.top:
#include "IBU.acpype/atomtypes.itp"
#include "protein.itp"
#include "IBU.acpype/IBU_GMX.itp"

[ molecules ]
Protein_chain_A    1
IBU                1

# 5. Solvate and continue normal workflow
gmx_mpi editconf -f complex.gro -o boxed.gro -c -d 1.2 -bt dodecahedron
gmx_mpi solvate -cp boxed.gro -o solvated.gro -p topol.top
```

## Advanced: Custom Parameters

For specialized molecules:

1. **Quantum chemistry:** Calculate charges (RESP, AM1-BCC)
2. **Bonded parameters:** Derive from Hessian (vibrational frequencies)
3. **Validation:** Compare with QM energies, experimental data

See force field papers for parameter derivation protocols.

## Resources

- **ACPYPE:** https://github.com/alanwilter/acpype
- **CGenFF:** https://cgenff.umaryland.edu/
- **ATB:** https://atb.uq.edu.au/
- **LigParGen:** http://zarbi.chem.yale.edu/ligpargen/
- **RDKit:** https://www.rdkit.org/docs/
- **PubChem:** https://pubchem.ncbi.nlm.nih.gov/
