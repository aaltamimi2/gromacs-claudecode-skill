# ORCA COSMO-RS Workflow Guide

Complete toolkit for converting Gaussian workflows to ORCA and running COSMO-RS calculations.

## Quick Start

```bash
# Method 1: From Gaussian .gjf file
python gjf_to_orca.py your_file.gjf
./orca_cosmo_workflow.sh your_file.xyz 56

# Method 2: From PDB file
python pdb_to_xyz.py molecule.pdb
./orca_cosmo_workflow.sh molecule.xyz 56

# Method 3: Direct from XYZ file
./orca_cosmo_workflow.sh molecule.xyz 56
```

## Workflow Comparison

### Original Gaussian Workflow
```
1. # bvp86/tzvp/dga1 scf=tight geom=connectivity opt=tight
   → Gas phase optimization

2. # bvp86/tzvp/dga1 opt=tight geom=checkpoint scf=tight SCRF=(CPCM)
   → CPCM optimization in Water

3. # bvp86/tzvp/dga1 geom=checkpoint guess=read scf=tight SCRF=COSMORS
   → COSMO-RS single point (infinite dielectric)
```

### Equivalent ORCA Workflow
```
1. ! BP86 def2-TZVP def2/J TightSCF TightOpt Grid5 RI-J
   → Gas phase optimization

2. ! BP86 def2-TZVP def2/J TightSCF TightOpt Grid5 CPCM(Water) RI-J
   → CPCM optimization in Water

3. ! BP86 def2-TZVP def2/J TightSCF Grid5 RI-J
   %cpcm epsilon infinity end
   → COSMO-RS single point (infinite dielectric)
```

## Scripts Overview

### 1. gjf_to_orca.py
Converts Gaussian .gjf files to ORCA input format.

**Usage:**
```bash
python gjf_to_orca.py input.gjf [output_name]
```

**Features:**
- Parses Gaussian route section
- Handles PDB-style atom names
- Extracts charge and multiplicity
- Creates both .inp and .xyz files
- Preserves large coordinates (for polymers)

**Example:**
```bash
python gjf_to_orca.py PU_conf_2300.gjf conf_2300
# Creates: conf_2300.inp and conf_2300.xyz
```

### 2. pdb_to_xyz.py
Converts PDB files to XYZ format for ORCA.

**Usage:**
```bash
python pdb_to_xyz.py molecule.pdb [-c CHARGE] [-m MULT] [-o OUTPUT]
```

**Options:**
- `-c, --charge`: Molecular charge (default: 0)
- `-m, --mult`: Spin multiplicity (default: 1)
- `-o, --output`: Output XYZ filename
- `-t, --title`: Custom title line

**Examples:**
```bash
# Neutral molecule
python pdb_to_xyz.py molecule.pdb

# Cation, doublet
python pdb_to_xyz.py molecule.pdb -c 1 -m 2

# Custom output
python pdb_to_xyz.py molecule.pdb -o my_structure.xyz
```

### 3. orca_cosmo_workflow.sh
Complete 3-step COSMO-RS workflow (gas → CPCM → COSMO).

**Usage:**
```bash
./orca_cosmo_workflow.sh molecule.xyz [nprocs]
```

**Parameters:**
- `molecule.xyz`: Input XYZ file
- `nprocs`: Number of cores (default: 56)

**What it does:**
1. **Gas phase optimization** (BP86/def2-TZVP)
   - Creates: `molecule_gas.inp`, `molecule_gas.out`, `molecule_gas.xyz`
   
2. **CPCM optimization** in water (BP86/def2-TZVP)
   - Uses gas phase geometry as starting point
   - Creates: `molecule_cpcm.inp`, `molecule_cpcm.out`, `molecule_cpcm.xyz`
   
3. **COSMO-RS single point** (infinite dielectric)
   - Uses CPCM geometry
   - Creates: `molecule_cosmo.inp`, `molecule_cosmo.out`
   - **Generates: `molecule_cosmo.cosmo`** ← IMPORTANT FILE

**Output:**
- All input/output files for each step
- `.cosmo` file for use with COSMOtherm
- Energy summary and solvation free energy

## Complete Example Workflow

### Starting from Gaussian .gjf file:

```bash
# 1. Convert Gaussian file
python gjf_to_orca.py PU_conf_2300.gjf
# Creates: PU_conf_2300.inp, PU_conf_2300.xyz

# 2. Run COSMO-RS workflow
./orca_cosmo_workflow.sh PU_conf_2300.xyz 56

# 3. Wait for completion (check progress)
tail -f PU_conf_2300_gas.out      # Monitor gas phase
tail -f PU_conf_2300_cpcm.out     # Monitor CPCM
tail -f PU_conf_2300_cosmo.out    # Monitor COSMO

# 4. Check results
ls -lh PU_conf_2300_cosmo.cosmo   # Your COSMO file!
```

### Starting from PDB file:

```bash
# 1. Convert PDB to XYZ
python pdb_to_xyz.py molecule.pdb -c 0 -m 1

# 2. Run COSMO-RS workflow
./orca_cosmo_workflow.sh molecule.xyz 56

# 3. Get COSMO file
ls -lh molecule_cosmo.cosmo
```

## Batch Processing Multiple Molecules

```bash
#!/bin/bash
# Process multiple Gaussian files

for gjf in *.gjf; do
    name=$(basename "$gjf" .gjf)
    echo "Processing $name..."
    
    # Convert
    python gjf_to_orca.py "$gjf" "$name"
    
    # Run workflow
    ./orca_cosmo_workflow.sh "${name}.xyz" 56
    
    echo "Completed $name"
    echo "---"
done

# Collect all COSMO files
mkdir -p cosmo_files
cp *_cosmo.cosmo cosmo_files/
```

## File Format Details

### XYZ File Format
```
155
Polyurethane oligomer, charge=0 mult=1
O    127.51000000   18.24000000   42.61000000
C    127.97000000   17.55000000   41.46000000
...
```
- Line 1: Number of atoms
- Line 2: Comment (includes charge and multiplicity)
- Lines 3+: Element X Y Z

### COSMO File
The `.cosmo` file is a text file containing:
- Molecular geometry
- Segment information for COSMO surface
- Screening charges
- Surface areas

**Used by:** COSMOtherm, COSMOquick for property predictions

## Performance Notes

For a typical polyurethane oligomer (~155 atoms):
- **Gas phase opt**: 30-60 minutes (56 cores)
- **CPCM opt**: 45-90 minutes (56 cores)
- **COSMO SP**: 10-20 minutes (56 cores)
- **Total**: ~1.5-3 hours per molecule

## Troubleshooting

### "ORCA TERMINATED ABNORMALLY"

**Check:**
1. SCF convergence issues
   - Solution: Add `SlowConv` keyword or increase MaxIter
   
2. Optimization not converging
   - Solution: Use looser thresholds or start from better geometry
   
3. Memory issues
   - Solution: Reduce `%maxcore` or increase available RAM

### COSMO file not generated

**Check:**
1. COSMO-RS calculation completed successfully?
   ```bash
   grep "ORCA TERMINATED NORMALLY" molecule_cosmo.out
   ```

2. `%cpcm epsilon infinity` block present in input?
   ```bash
   grep -A 2 "%cpcm" molecule_cosmo.inp
   ```

### Large coordinate values (127-138 Å)

This is normal for structures from PDB files. ORCA handles absolute coordinates correctly.

## Comparison: ORCA vs Gaussian

| Feature | ORCA | Gaussian |
|---------|------|----------|
| License | Free (academic) | Commercial |
| Command | Single: `orca input.inp` | Complex setup |
| XYZ I/O | Native | Need conversion |
| COSMO files | Automatic | Automatic |
| Scripting | Very easy | Complex |
| Speed | Similar (RI methods) | Similar |
| Parallelization | Excellent | Good |

## Next Steps

After generating `.cosmo` files:

1. **COSMOtherm**: Use for property predictions
   - Solubility
   - Partition coefficients
   - Activity coefficients
   - VLE data

2. **Machine Learning**: Use COSMO-RS descriptors
   - σ-profiles
   - σ-moments
   - Surface areas

## Additional Resources

- ORCA Manual: https://www.orcasoftware.de/tutorials_orca/
- ORCA Forum: https://orcaforum.kofo.mpg.de
- COSMO-RS Theory: COSMOtherm documentation

## Notes

- BP86 ≈ BVP86 (same functional)
- def2-TZVP ≈ TZVP (same basis set, better organized)
- RI-J provides ~2-3× speedup with minimal accuracy loss
- All scripts preserve exact methodology from Gaussian workflow
