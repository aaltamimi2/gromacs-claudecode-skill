# GROMACS Troubleshooting Guide

Common errors, their causes, and solutions.

## LINCS Warnings and Errors

### Symptom
```
LINCS WARNING: relative constraint deviation after LINCS: 0.5000
```
Or fatal: `Too many LINCS warnings`

### Causes
1. Bad initial contacts (overlapping atoms)
2. Forces too large for timestep
3. Inappropriate constraints for the system
4. Poor equilibration

### Solutions

**Immediate fixes:**
```bash
# Run longer/stricter energy minimization
gmx grompp -f em_strict.mdp -c system.gro -p topol.top -o em.tpr
# em_strict.mdp: emtol = 100, nsteps = 100000
```

**Reduce timestep:**
```mdp
dt = 0.001   ; 1 fs instead of 2 fs
```

**Increase LINCS iterations:**
```mdp
lincs_iter = 2   ; Default is 1
lincs_order = 6  ; Default is 4
```

**Check for bad contacts:**
```bash
gmx editconf -f system.gro -o system.pdb
# Visual inspection in VMD/PyMOL
# Look for overlapping atoms, especially at interfaces
```

## System Explodes / Coordinates Blow Up

### Symptom
```
Fatal error: Atom X in domain Y moved more than Z nm
```
Or trajectory shows atoms flying apart.

### Causes
1. Overlapping atoms in initial structure
2. Missing or incorrect topology
3. Wrong units (nm vs Å)
4. Incompatible force field / coordinate combination

### Solutions

**Check topology:**
```bash
gmx grompp -f em.mdp -c system.gro -p topol.top -o em.tpr -pp processed.top
# Examine processed.top for correct atom counts, missing parameters
```

**Verify coordinates:**
```bash
gmx editconf -f system.gro -o check.pdb
# Should be nm, typical protein ~5-10 nm box
```

**Visual inspection:**
```bash
vmd system.gro  # or PyMOL
# Check: no overlapping atoms, reasonable box size, correct periodicity
```

**Gradual equilibration:**
```bash
# Shorter timestep initially
dt = 0.0005  ; 0.5 fs for initial equilibration
# Then gradually increase to 2 fs
```

## NaN in Energy / Forces

### Symptom
```
Step X: NaN in energy
Fatal error: NaN or Inf in forces
```

### Causes
1. Numerical instability (often from LINCS issues)
2. Division by zero (atoms too close)
3. Corrupted checkpoint
4. PME grid issues

### Solutions

**Start fresh from last good checkpoint:**
```bash
gmx check -f state.cpt  # Verify checkpoint integrity
gmx mdrun -cpi state_previous.cpt -deffnm md
```

**Check PME settings:**
```mdp
fourierspacing = 0.12   ; Try 0.16 for troubleshooting
pme_order = 4           ; Standard value
```

**Reduce timestep and re-equilibrate:**
```bash
# Go back to NPT with smaller dt
dt = 0.001
nsteps = 100000
```

## Density Issues (Too Low/High)

### Symptom
- Density far from expected (~1000 kg/m³ for water)
- Box shrinks/expands excessively

### Causes
1. Wrong pressure coupling parameters
2. Incompatible compressibility
3. Incorrect reference pressure
4. System not equilibrated

### Solutions

**Check barostat settings:**
```mdp
pcoupl          = C-rescale
pcoupltype      = isotropic
tau_p           = 2.0           ; Not too fast
ref_p           = 1.0           ; bar
compressibility = 4.5e-5        ; Water at 300K
```

**Verify density calculation:**
```bash
gmx energy -f npt.edr -o density.xvg
# Select "Density"
# Should stabilize around 1000 kg/m³ for aqueous
```

**Longer NPT equilibration:**
```bash
nsteps = 500000  ; 1 ns NPT before production
```

## Temperature Drift

### Symptom
- Temperature drifts from target
- Different groups at different temperatures

### Causes
1. tc-grps don't match index groups
2. Wrong tau_t coupling time
3. Flying ice cube effect (rare with v-rescale)

### Solutions

**Verify temperature groups:**
```bash
gmx make_ndx -f system.gro -o index.ndx
# Check available groups match tc-grps in .mdp
```

**Correct .mdp:**
```mdp
tcoupl      = v-rescale
tc-grps     = Protein Non-Protein   ; Must exist in index
tau_t       = 0.1     0.1
ref_t       = 300     300
```

**Check temperature output:**
```bash
gmx energy -f md.edr -o temp.xvg
# Select temperature terms for each group
```

## Pressure Oscillations

### Symptom
- Large pressure fluctuations
- Box size oscillates significantly

### Causes
1. Barostat coupling too tight
2. Small system size
3. Incompatible with thermostat

### Solutions

**Increase barostat coupling time:**
```mdp
tau_p = 5.0   ; Slower coupling (default often 2.0)
```

**Note:** Pressure fluctuations are normal. For small systems, instantaneous pressure can vary by hundreds of bar. Check the *running average*, not instantaneous values.

## grompp Errors

### "No default X type" / Missing parameters
```
Fatal error: No default Proper Dih. types
```

**Solution:** Force field missing parameters for your molecule.
```bash
# Check if residue is supported
grep "YOURRES" charmm36.ff/residuetypes.dat

# May need to generate topology separately
# Use CGenFF, ATB, or parameterize manually
```

### "Atom X in residue Y was not found"
**Solution:** Coordinate file atom names don't match topology.
```bash
# Compare:
grep "ATOM" system.pdb | head
grep "\[ atoms \]" -A 20 topol.top
```

### "System has non-zero total charge"
**Solution:** Add ions to neutralize.
```bash
gmx genion -s ions.tpr -o system.gro -p topol.top -pname NA -nname CL -neutral
```

## Performance Issues

### Slow simulation
1. **Check load balance:**
```bash
tail -50 md.log | grep -A 20 "Performance"
```

2. **GPU utilization:**
```bash
nvidia-smi  # Check GPU usage during run
```

3. **Optimize settings:**
```mdp
nstlist = 40        ; Increase for better GPU performance
```

4. **MPI/OpenMP balance:**
```bash
# Typical: 1 MPI rank per GPU, fill OMP threads
gmx mdrun -ntmpi 1 -ntomp 16 -nb gpu -pme gpu -bonded gpu
```

### "Domain decomposition failed"
**Solution:** Try different decomposition or single rank:
```bash
gmx mdrun -dd 2 2 2 ...  # Explicit decomposition
# Or for small systems:
gmx mdrun -ntmpi 1 ...   # Single rank
```

## Checkpoint/Restart Issues

### "Checkpoint mismatch"
**Solution:** Checkpoint and .tpr must match exactly.
```bash
# Re-run grompp with same parameters
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md -cpi state.cpt
```

### Extending simulations
```bash
# Method 1: convert-tpr (preferred for simple extension)
gmx convert-tpr -s md.tpr -extend 10000 -o md_extended.tpr
gmx mdrun -s md_extended.tpr -cpi state.cpt -deffnm md -noappend

# Method 2: New grompp (if changing parameters)
gmx grompp -f md_new.mdp -c md.gro -t state.cpt -p topol.top -o md_new.tpr
```

## Quick Diagnostic Commands

```bash
# Check trajectory integrity
gmx check -f traj.xtc

# Examine energy file
gmx energy -f ener.edr  # Interactive selection

# Check .tpr contents
gmx dump -s run.tpr | head -200

# Compare .tpr files
gmx check -s1 old.tpr -s2 new.tpr

# Trajectory info
gmx check -f traj.xtc 2>&1 | grep -E "^(Coords|Step|Time)"
```
