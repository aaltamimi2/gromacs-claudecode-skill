# GROMACS MPI Configuration

## Understanding gmx vs gmx_mpi

GROMACS can be compiled in two modes:

1. **gmx** - Thread-MPI version (easier for single-node, no mpirun required)
2. **gmx_mpi** - Full MPI version (requires mpirun/mpiexec wrapper)

## Detecting Your Installation

```bash
# Check which version you have
which gmx       # Thread-MPI version
which gmx_mpi   # Full MPI version

# Check version details
gmx --version 2>&1 | grep -i mpi
gmx_mpi --version 2>&1 | head -3
```

## Running Commands with gmx_mpi

If you have `gmx_mpi`, you **must** use `mpirun` or `mpiexec`:

### Single Process (Development/Testing)
```bash
mpirun -np 1 gmx_mpi pdb2gmx -f input.pdb -o output.gro
mpirun -np 1 gmx_mpi mdrun -deffnm md
```

### Multi-Process (Production)
```bash
# CPU cluster (16 MPI ranks)
mpirun -np 16 gmx_mpi mdrun -deffnm md

# GPU cluster (1 rank per GPU, fill OpenMP threads)
mpirun -np 4 gmx_mpi mdrun -deffnm md -ntomp 8 -nb gpu -pme gpu
```

## Common MPI Errors

### Error: "write_line error; fd=-1"
```
application called MPI_Abort(MPI_COMM_WORLD, 1) - process 0
[unset]: write_line error; fd=-1 buf=:cmd=abort exitcode=1
```

**Cause**: Running `gmx_mpi` directly without mpirun

**Fix**: Add `mpirun -np 1` prefix
```bash
# Wrong:
gmx_mpi pdb2gmx -f input.pdb

# Correct:
mpirun -np 1 gmx_mpi pdb2gmx -f input.pdb
```

### Error: "MPI_Init not called"

**Fix**: Same as above - use mpirun wrapper

## Using the gmx_wrapper.sh Script

The skill includes a wrapper that automatically handles this:

```bash
# Automatically detects gmx vs gmx_mpi
./scripts/gmx_wrapper.sh pdb2gmx -f input.pdb -o output.gro
./scripts/gmx_wrapper.sh mdrun -deffnm md
```

To use it system-wide:
```bash
# Add alias to your shell profile
alias gmx='/path/to/gromacs-skill/scripts/gmx_wrapper.sh'
```

## HPC Environments

On clusters, GROMACS is often in modules:

```bash
# Load GROMACS module
module load gromacs/2021.7
module load gromacs/2024.1-gpu

# Check what's available
module avail gromacs

# Some clusters provide both versions
which gmx       # Thread-MPI
which gmx_mpi   # Full MPI
```

### SLURM Job Example with gmx_mpi

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00

module load gromacs/2021.7

# Use srun (SLURM's mpirun equivalent)
srun gmx_mpi mdrun -deffnm production
```

## Performance Recommendations

### Thread-MPI (gmx)
```bash
# Let GROMACS auto-detect optimal setup
gmx mdrun -deffnm md

# Or specify explicitly
gmx mdrun -deffnm md -ntmpi 4 -ntomp 8 -nb gpu
```

### Full MPI (gmx_mpi)
```bash
# CPU-only: 1 rank per core
mpirun -np 32 gmx_mpi mdrun -deffnm md -ntomp 1

# GPU: 1 rank per GPU, fill OpenMP
mpirun -np 4 gmx_mpi mdrun -deffnm md -ntomp 16 -nb gpu -pme gpu
```

## Troubleshooting

### Check MPI is working
```bash
mpirun -np 2 hostname
# Should print hostname twice
```

### Debug MPI issues
```bash
# Verbose MPI output
mpirun -np 1 -v gmx_mpi --version

# Check MPI library
ldd $(which gmx_mpi) | grep mpi
```

### Common solutions
1. **Module not loaded**: `module load openmpi` or `module load gromacs`
2. **Wrong MPI library**: Recompile GROMACS with correct MPI
3. **Firewall blocking**: Check MPI communication ports
4. **Hybrid builds**: Some installations have both `gmx` and `gmx_mpi` - prefer `gmx` for simplicity
