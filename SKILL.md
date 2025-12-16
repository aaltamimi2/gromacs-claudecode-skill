---
name: gromacs
description: "Create, run, and analyze molecular dynamics simulations using GROMACS. Use this skill when the user asks for: (1) System setup and topology preparation, (2) .mdp parameter file design, (3) Simulation workflows (EM → NVT → NPT → production), (4) HPC job scripts and performance tuning, (5) Trajectory analysis (energy, density, MSD, SASA, RDF), (6) Free energy calculations (umbrella sampling, FEP, PMF), (7) Troubleshooting simulation failures (LINCS, exploding systems, barostat instability)."
---

# GROMACS Molecular Dynamics Skill

Guide for reproducible, production-grade GROMACS MD workflows: system building → equilibration → production → analysis.

## Design Thinking

Before writing commands, identify:
- **Goal**: Stability, sampling, kinetics, thermodynamics, free energy?
- **Model**: Force field, water model, constraints strategy
- **Risk**: LINCS failures, bad contacts, barostat instability, PBC artifacts
- **Deliverable**: Command sequence, job script, .mdp blocks, analysis recipe

## Standard Workflow

### 1. System Preparation
```bash
# Generate topology (example for protein)
gmx pdb2gmx -f input.pdb -o processed.gro -water tip3p -ff charmm27

# Define box and solvate
gmx editconf -f processed.gro -o boxed.gro -c -d 1.2 -bt dodecahedron
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

# Add ions to neutralize
gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o system.gro -p topol.top -pname NA -nname CL -neutral
```

### 2. Energy Minimization
```bash
gmx grompp -f em.mdp -c system.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
```
**Validation**: Check `em.log` for `Fmax < emtol` and `potential energy < 0`.

### 3. Equilibration (NVT → NPT)
```bash
# NVT (temperature coupling)
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

# NPT (pressure coupling)
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt
```
**Validation**: Temperature stable at target ±5K, density ~1000 kg/m³ for aqueous.

### 4. Production MD
```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md
```

## Core Tools Reference

### gmx grompp (Preprocessor)
Creates `.tpr` run input from topology + coordinates + parameters.

Key flags:
- `-f` .mdp parameters, `-p` topology, `-c` coordinates
- `-t` checkpoint for continuation (velocities, coupling variables)
- `-r` restraint reference coordinates (required for position restraints)
- `-n` custom index groups
- `-po mdout.mdp` shows actual interpreted parameters
- `-pp processed.top` shows expanded topology (debug macros/includes)

### gmx mdrun (Engine)
Executes simulation from `.tpr`.

Key flags:
- `-deffnm` sets base name for all outputs
- `-cpi state.cpt` continues from checkpoint
- `-maxh N` wall-time limit (writes checkpoint before exit)
- `-nb gpu -pme gpu -bonded gpu -update gpu` GPU acceleration
- `-ntmpi N -ntomp M` MPI ranks × OpenMP threads

### Analysis Tools

| Tool | Purpose | Key Options |
|------|---------|-------------|
| `gmx energy` | Energy/pressure/temperature extraction | Interactive selection |
| `gmx density` | Density profiles along axis | `-center`, `-symm` for bilayers |
| `gmx msd` | Mean-square displacement, diffusion | `-mol` for per-molecule, `-beginfit/-endfit` |
| `gmx sasa` | Solvent-accessible surface area | `-surface`, `-output` selections |
| `gmx rdf` | Radial distribution functions | `-ref`, `-sel` groups |
| `gmx trjconv` | Trajectory processing/PBC fixing | `-pbc mol`, `-center`, `-fit` |

## Default Parameters

Unless specified otherwise, use:
- **Thermostat**: `v-rescale` (stochastic velocity rescaling, tau_t = 0.1-1.0 ps)
- **Barostat**: `C-rescale` (stochastic cell rescaling, tau_p = 2.0-5.0 ps)
- **Cutoffs**: 1.2 nm for vdW, PME for electrostatics
- **Constraints**: `h-bonds` with LINCS (or `all-bonds` if needed)
- **Time step**: 2 fs (with h-bond constraints), 1 fs (unconstrained)

## Reference Files

For detailed information, see:
- **[references/mdp-templates.md](references/mdp-templates.md)**: Complete .mdp examples for EM, NVT, NPT, production, and specialized runs
- **[references/free-energy.md](references/free-energy.md)**: Umbrella sampling, steered MD, PMF calculations, FEP
- **[references/troubleshooting.md](references/troubleshooting.md)**: Common errors (LINCS, exploding systems, NaN) and fixes
- **[references/hpc-scripts.md](references/hpc-scripts.md)**: SLURM job templates for CPU/GPU clusters
- **[references/mpi-configuration.md](references/mpi-configuration.md)**: gmx vs gmx_mpi setup, PLUMED issues

## Quick Troubleshooting

| Symptom | Likely Cause | Fix |
|---------|--------------|-----|
| LINCS warning | Bad contacts, large forces | Longer/stricter EM, reduce dt |
| System explodes | Overlapping atoms, bad topology | Check topology, visual inspection |
| NaN in energy | Unstable simulation | Reduce dt, check parameters |
| Density too low/high | Wrong pressure coupling | Check barostat settings, ref_p |
| Temperature drift | Thermostat misconfigured | Check tc-grps match index |
| PLUMED symbol error | PLUMED library not loaded | Set LD_LIBRARY_PATH or use gmx |

See [references/troubleshooting.md](references/troubleshooting.md) and [references/mpi-configuration.md](references/mpi-configuration.md) for detailed diagnostics.

**Note**: Examples use `gmx` for simplicity. Prefer `gmx_mpi` when available (see mpi-configuration.md).
