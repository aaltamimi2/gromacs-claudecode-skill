# GROMACS Claude Code Skill

A comprehensive Claude Code skill for molecular dynamics simulations using GROMACS.

## Installation

Clone this repository and the skill will be automatically detected by Claude Code from `.claude/skills/gromacs-skill/`:

```bash
git clone https://github.com/aaltamimi2/gromacs-claudecode-skill.git
cd gromacs-claudecode-skill
```

The `.claude/skills/gromacs-skill/` directory contains symlinks to the root-level SKILL.md, references/, and scripts/ directories.

## Usage

Once registered, you can invoke the skill by asking Claude for GROMACS-related tasks:

- "Help me set up a protein MD simulation"
- "Create an energy minimization workflow"
- "Design .mdp parameters for NPT equilibration"
- "Analyze trajectory files for density and temperature"
- "Debug LINCS warnings in my simulation"
- "Set up umbrella sampling for free energy calculations"

## Features

- **System Preparation**: Topology generation, solvation, ion addition
- **Simulation Workflows**: EM → NVT → NPT → Production pipelines
- **Parameter Files**: Comprehensive .mdp templates for all simulation types
- **HPC Integration**: SLURM job scripts, GPU optimization
- **Analysis Tools**: Energy, density, RDF, SASA, MSD calculations
- **Free Energy**: Umbrella sampling, FEP, PMF calculations
- **Troubleshooting**: Solutions for LINCS, exploding systems, NaN errors, MPI issues

## GROMACS Configuration

This skill supports both `gmx_mpi` (preferred) and `gmx` (fallback).

**Workstation usage** - Run directly without mpirun:
```bash
gmx_mpi pdb2gmx -f input.pdb -o output.gro
gmx_mpi mdrun -deffnm production
```

**HPC cluster usage** - Use job scheduler launcher:
```bash
srun gmx_mpi mdrun -deffnm production
# or: mpirun -np 16 gmx_mpi mdrun -deffnm production
```

See `references/mpi-configuration.md` for detailed configuration and PLUMED troubleshooting.

### Using the Wrapper Script

A wrapper script is provided to automatically handle MPI:

```bash
scripts/gmx_wrapper.sh pdb2gmx -f input.pdb
```

Or create an alias:
```bash
alias gmx='<full-path>/gromacs-claudecode-skill/scripts/gmx_wrapper.sh'
```

## Repository Structure

```
gromacs-claudecode-skill/
├── SKILL.md                          # Main skill file
├── references/                       # Detailed documentation
│   ├── mdp-templates.md             # .mdp parameter examples
│   ├── free-energy.md               # Free energy calculations
│   ├── troubleshooting.md           # Error diagnosis and fixes
│   ├── hpc-scripts.md               # HPC job templates
│   └── mpi-configuration.md         # MPI setup guide
├── scripts/                          # Utility scripts
│   ├── gmx_wrapper.sh               # Automatic gmx/gmx_mpi handler
│   └── check_equilibration.py       # Equilibration validation
└── .claude/skills/gromacs-skill/    # Symlinks for Claude Code
    ├── SKILL.md -> ../../../SKILL.md
    ├── references -> ../../../references
    └── scripts -> ../../../scripts
```

## Quick Start Example

```bash
# 1. Generate topology
gmx_mpi pdb2gmx -f protein.pdb -o processed.gro -ff amber99sb-ildn

# 2. Create box and solvate
gmx_mpi editconf -f processed.gro -o boxed.gro -c -d 1.2
gmx_mpi solvate -cp boxed.gro -o solvated.gro -p topol.top

# 3. Energy minimization
gmx_mpi grompp -f em.mdp -c solvated.gro -p topol.top -o em.tpr
gmx_mpi mdrun -deffnm em

# 4. NVT equilibration
gmx_mpi grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
gmx_mpi mdrun -deffnm nvt

# 5. Check equilibration
python scripts/check_equilibration.py nvt.edr
```

## Troubleshooting Common Errors

### MPI write_line error
**Error**: `write_line error; fd=-1 buf=:cmd=abort exitcode=1`

**Solution**: Rare MPI configuration issue - try using `gmx` instead of `gmx_mpi`

### PLUMED symbol error
**Error**: `gmx_mpi: symbol lookup error: undefined symbol: plumed_hrex`

**Solution**: Load PLUMED library - `export LD_LIBRARY_PATH=/path/to/plumed/lib:$LD_LIBRARY_PATH`

### LINCS warnings
**Error**: `LINCS WARNING relative constraint deviation`

**Solution**: Longer energy minimization, reduce timestep to 1 fs

See the troubleshooting reference for comprehensive error solutions.

## Contributing

To extend the skill:
1. Edit `SKILL.md` for main instructions
2. Add detailed guides in `references/`
3. Add utility scripts in `scripts/`

The `.claude/skills/gromacs-skill/` directory will automatically pick up changes via symlinks.

## License

This skill provides guidance for using GROMACS, which is licensed under LGPL.
