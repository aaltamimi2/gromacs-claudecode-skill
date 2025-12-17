# GROMACS Configuration Guide

## Important: Run Commands As-Is

**If your environment is properly configured, run `gmx_mpi` commands exactly as specified:**

```bash
# ✓ Correct - run directly
gmx_mpi insert-molecules -f input.gro -ci ligand.gro -nmol 100 -o output.gro
gmx_mpi grompp -f md.mdp -c input.gro -p topol.top -o md.tpr
gmx_mpi mdrun -deffnm md
```

**DO NOT add workarounds unless specifically needed:**
```bash
# ✗ Wrong - unnecessary modifications
unset PLUMED_KERNEL && gmx_mpi ...    # Don't add this
export LD_LIBRARY_PATH=... && gmx_mpi ... # Don't add this inline
```

## Environment Setup for Claude Code

**Claude Code runs bash in non-interactive mode**, which doesn't source `~/.bashrc` automatically. Set environment variables system-wide:

### Option 1: Use ~/.bash_profile (Recommended)

```bash
# Add to ~/.bash_profile (loaded for all sessions)
export GMXLIB=/home/aaltamimi2/gromacs-plumed/share/gromacs/top
export LD_LIBRARY_PATH=/home/aaltamimi2/gromacs-plumed/lib:$LD_LIBRARY_PATH
export PATH=/home/aaltamimi2/gromacs-plumed/bin:$PATH
```

### Option 2: Use /etc/environment (System-wide)

```bash
# Add to /etc/environment (requires sudo)
GMXLIB=/home/aaltamimi2/gromacs-plumed/share/gromacs/top
LD_LIBRARY_PATH=/home/aaltamimi2/gromacs-plumed/lib:$LD_LIBRARY_PATH
PATH=/home/aaltamimi2/gromacs-plumed/bin:$PATH
```

### Option 3: Source in ~/.bashrc AND export for subshells

```bash
# Add to ~/.bashrc
export GMXLIB=/home/aaltamimi2/gromacs-plumed/share/gromacs/top
export LD_LIBRARY_PATH=/home/aaltamimi2/gromacs-plumed/lib:$LD_LIBRARY_PATH
export PATH=/home/aaltamimi2/gromacs-plumed/bin:$PATH

# Also source ~/.bashrc in ~/.bash_profile
echo 'source ~/.bashrc' >> ~/.bash_profile
```

**After setup, restart your shell or run:**
```bash
source ~/.bash_profile
```

Once configured, commands work directly without prefixes.

## Understanding gmx vs gmx_mpi

GROMACS can be compiled in two modes:

1. **gmx_mpi** - Full MPI version (preferred for flexibility)
2. **gmx** - Thread-MPI version (fallback option)

**Recommended**: Use `gmx_mpi` when available, fall back to `gmx` if not.

## Detecting Your Installation

```bash
# Check which versions you have (priority order)
which gmx_mpi   # Preferred
which gmx       # Fallback

# Check version details
gmx_mpi --version 2>&1 | head -3
gmx --version 2>&1 | head -3
```

## Running Commands

### Workstation Usage (Single Node)

On a workstation, you can run `gmx_mpi` **directly without mpirun**:

```bash
# All commands work directly
gmx_mpi pdb2gmx -f input.pdb -o output.gro
gmx_mpi editconf -f input.gro -o boxed.gro -box 5 5 5
gmx_mpi grompp -f md.mdp -c input.gro -p topol.top -o md.tpr
gmx_mpi mdrun -deffnm md
```

GROMACS will automatically use thread-MPI mode when called this way.

### HPC/Cluster Usage (Multi-Node)

On HPC systems, use the job scheduler's MPI launcher:

```bash
# SLURM
srun gmx_mpi mdrun -deffnm production

# Or with explicit MPI
mpirun -np 16 gmx_mpi mdrun -deffnm production
```

## Common Errors and Solutions

### Error: PLUMED Symbol Lookup Error

```
gmx_mpi: symbol lookup error: gmx_mpi: undefined symbol: plumed_hrex
```

**Cause**: Your `gmx_mpi` was compiled with PLUMED support but the PLUMED library isn't loaded.

**Solution 1: Load PLUMED Module** (HPC)
```bash
module load plumed
# Or
module load plumed/2.9

# Verify
ldd $(which gmx_mpi) | grep plumed
```

**Solution 2: Set LD_LIBRARY_PATH** (Workstation - One-Time Setup)

**Add to `~/.bashrc` for permanent fix:**
```bash
# Add PLUMED library path permanently
echo 'export LD_LIBRARY_PATH=/path/to/plumed/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc
```

**Example:**
```bash
# Find PLUMED library first
find /home/$(whoami) -name "libplumedKernel.so*" 2>/dev/null

# Then add to ~/.bashrc
echo 'export LD_LIBRARY_PATH=/home/username/gromacs-plumed/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc
```

**DO NOT use inline environment modifications:**
```bash
# ✗ Wrong - adds complexity to every command
export LD_LIBRARY_PATH=/path/to/lib:$LD_LIBRARY_PATH && gmx_mpi ...
unset PLUMED_KERNEL && gmx_mpi ...

# ✓ Correct - set once in ~/.bashrc, then use commands normally
gmx_mpi ...
```

**Solution 3: Find PLUMED Installation**
```bash
# Locate PLUMED
find /usr -name "libplumedKernel.so*" 2>/dev/null
find /opt -name "libplumedKernel.so*" 2>/dev/null
ldconfig -p | grep plumed

# Once found, add to LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/path/to/plumed/lib:$LD_LIBRARY_PATH
```

**Solution 4: Use Non-PLUMED GROMACS**
```bash
# If you don't need PLUMED, use gmx instead
which gmx
gmx --version  # Check if available

# Or install/compile GROMACS without PLUMED
```

### Error: "write_line error; fd=-1"

```
application called MPI_Abort(MPI_COMM_WORLD, 1) - process 0
[unset]: write_line error; fd=-1 buf=:cmd=abort exitcode=1
```

**Cause**: Old MPI configuration issue (rare with modern GROMACS)

**Fix**: This typically indicates a system MPI issue. Try:
```bash
# Check MPI installation
which mpirun
mpirun --version

# Test MPI
mpirun -np 2 hostname

# If MPI is broken, use gmx (thread-MPI) instead of gmx_mpi
```

## Wrapper Script

The repository includes `scripts/gmx_wrapper.sh` that automatically selects `gmx_mpi` (preferred) or `gmx` (fallback):

```bash
# Use wrapper
~/gromacs-claudecode-skill/scripts/gmx_wrapper.sh pdb2gmx -f input.pdb

# Or create alias
alias gmx='~/gromacs-claudecode-skill/scripts/gmx_wrapper.sh'
```

The wrapper prioritizes `gmx_mpi` without requiring mpirun on workstations.

## HPC Environments

On clusters, load GROMACS via modules:

```bash
# Load GROMACS
module load gromacs

# Or specific version
module load gromacs/2023.1
module load gromacs/2021.7-plumed  # With PLUMED

# Check available versions
module avail gromacs

# Verify what's loaded
module list
which gmx_mpi
```

### SLURM Job Script Example

```bash
#!/bin/bash
#SBATCH --job-name=md_production
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --partition=compute

# Load modules
module load gromacs/2023.1

# Run with SLURM's srun (preferred on SLURM)
srun gmx_mpi mdrun -deffnm production

# Or with mpirun
# mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm production
```

### GPU Job Script

```bash
#!/bin/bash
#SBATCH --job-name=md_gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --time=48:00:00

module load gromacs/2023.1-gpu
module load cuda/12.0

# 1 MPI rank per GPU, fill OpenMP threads
srun gmx_mpi mdrun -deffnm md -ntomp 16 -nb gpu -pme gpu -bonded gpu
```

## Performance Recommendations

### Workstation (Thread-MPI mode)
```bash
# Let GROMACS auto-detect
gmx_mpi mdrun -deffnm md

# Or specify threads
gmx_mpi mdrun -deffnm md -ntomp 16

# With GPU
gmx_mpi mdrun -deffnm md -ntomp 16 -nb gpu -pme gpu
```

### HPC Cluster (Full MPI mode)
```bash
# CPU-only: maximize MPI ranks
srun -n 128 gmx_mpi mdrun -deffnm md -ntomp 1

# GPU: 1 rank per GPU, fill OpenMP
srun -n 4 gmx_mpi mdrun -deffnm md -ntomp 32 -nb gpu -pme gpu
```

## Checking Your Installation

```bash
# 1. Which GROMACS do you have?
which gmx_mpi gmx

# 2. Check version and features
gmx_mpi --version

# 3. Check for PLUMED
gmx_mpi --version 2>&1 | grep -i plumed
ldd $(which gmx_mpi) | grep plumed

# 4. Check MPI library
ldd $(which gmx_mpi) | grep mpi

# 5. Test run
echo "Integrator = steep" > test.mdp
gmx_mpi grompp -f test.mdp -c /dev/null -p /dev/null -o test.tpr 2>&1 | head
```

## Setup for Your Environment

Add to `~/.bashrc`:

```bash
# GROMACS setup
export GMXLIB=/path/to/gromacs/share/gromacs/top

# For PLUMED-enabled GROMACS (if needed)
export LD_LIBRARY_PATH=/path/to/plumed/lib:$LD_LIBRARY_PATH

# Alias for convenience
alias gmx='gmx_mpi'  # Use gmx_mpi by default

# Or use wrapper
alias gmx='~/gromacs-claudecode-skill/scripts/gmx_wrapper.sh'
```

## Summary

| Environment | Command | Notes |
|-------------|---------|-------|
| **Workstation** | `gmx_mpi` | Direct call, no mpirun needed |
| **HPC (SLURM)** | `srun gmx_mpi` | Use job scheduler launcher |
| **HPC (Other)** | `mpirun -np N gmx_mpi` | Explicit MPI launcher |
| **Fallback** | `gmx` | Use if gmx_mpi unavailable |

**Priority**: gmx_mpi > gmx

**PLUMED**: If gmx_mpi compiled with PLUMED, ensure PLUMED library is in LD_LIBRARY_PATH
