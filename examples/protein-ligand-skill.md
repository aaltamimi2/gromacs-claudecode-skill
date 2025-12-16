---
name: protein-ligand-md
description: "Specialized skill for protein-ligand MD simulations. Covers ligand parameterization, complex preparation, binding site analysis, and binding free energy calculations."
---

# Protein-Ligand MD Simulation Skill

Extension of the GROMACS skill for protein-ligand systems.

## Ligand Parameterization

### CGenFF (CHARMM General Force Field)
```bash
# 1. Convert ligand to mol2
obabel ligand.pdb -O ligand.mol2

# 2. Use CGenFF server or local installation
# Upload to https://cgenff.umaryland.edu/
# Or use local:
cgenff ligand.mol2 > ligand.str

# 3. Convert to GROMACS format
python cgenff_charmm2gmx.py ligand ligand.mol2 ligand.str charmm36.ff
```

### ACPYPE (AMBER)
```bash
# Generate parameters using GAFF
acpype -i ligand.pdb -n 0  # 0 = net charge

# This creates:
# - ligand.acpype/ligand_GMX.gro
# - ligand.acpype/ligand_GMX.top
```

## Complex Preparation

### Combine Protein and Ligand
```bash
# 1. Process protein
gmx pdb2gmx -f protein.pdb -o protein.gro -p protein.top -ff charmm36

# 2. Combine structures
# Edit complex.gro manually or:
cat protein.gro > complex.gro
grep -v "Generated\|SOL\|#" ligand.gro >> complex.gro

# 3. Update topology
# Add ligand .itp to protein topology:
# In topol.top, add before [ system ]:
#include "ligand.itp"

# In [ molecules ] section:
# Protein_chain_A    1
# LIG                1
```

## Binding Site Analysis

### Key Interactions
```bash
# H-bonds between protein and ligand
gmx hbond -f md.xtc -s md.tpr -num hbnum.xvg

# Distance between specific atoms
gmx distance -f md.xtc -s md.tpr -select 'resname LIG and name C1' 'resname ASP and name OD1'

# Contacts
gmx mindist -f md.xtc -s md.tpr -od mindist.xvg -pi
```

## MM/PBSA or MM/GBSA

For binding free energy estimation:

```bash
# Use gmx_MMPBSA
gmx_MMPBSA -O -i mmpbsa.in -cs md.tpr -ci index.ndx -cg 1 13 -ct md.xtc
```

## Restraints for Complex Equilibration

```mdp
; Position restraints during initial equilibration
define = -DPOSRES -DPOSRES_LIG

; In topology, add restraint files:
; #ifdef POSRES_LIG
; #include "posre_ligand.itp"
; #endif
```

## Production MD Parameters

```mdp
; Extended sampling for binding modes
integrator      = md
dt              = 0.002
nsteps          = 250000000  ; 500 ns

; Enhanced output for analysis
nstxout-compressed = 5000     ; 10 ps frames
```

## Post-Analysis Workflow

1. **RMSD of protein and ligand separately**
2. **RMSF of binding site residues**
3. **Ligand orientation/pose stability**
4. **Protein-ligand interaction fingerprints**
5. **Binding free energy (MM/PBSA, MM/GBSA, or FEP)**
