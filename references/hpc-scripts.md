# HPC Job Script Templates

SLURM templates for running GROMACS on clusters.

## CPU-Only Job

```bash
#!/bin/bash
#SBATCH --job-name=gromacs_md
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Load modules (adjust for your cluster)
module load gromacs/2024.2

# Set OMP threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run simulation
gmx mdrun -deffnm md -ntomp $OMP_NUM_THREADS
```

## Single GPU Job

```bash
#!/bin/bash
#SBATCH --job-name=gromacs_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

module load gromacs/2024.2-cuda

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# GPU-accelerated run
gmx mdrun -deffnm md \
    -ntmpi 1 \
    -ntomp $OMP_NUM_THREADS \
    -nb gpu \
    -pme gpu \
    -bonded gpu \
    -update gpu
```

## Multi-GPU Job (Single Node)

```bash
#!/bin/bash
#SBATCH --job-name=gromacs_multi_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:4
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

module load gromacs/2024.2-cuda

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Multi-GPU: 1 MPI rank per GPU
gmx mdrun -deffnm md \
    -ntmpi 4 \
    -ntomp $OMP_NUM_THREADS \
    -nb gpu \
    -pme gpu \
    -npme 1
```

## Multi-Node MPI Job

```bash
#!/bin/bash
#SBATCH --job-name=gromacs_mpi
#SBATCH --partition=batch
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

module load gromacs/2024.2-mpi

# Large-scale MPI run
srun gmx_mpi mdrun -deffnm md
```

## Umbrella Sampling Array Job

```bash
#!/bin/bash
#SBATCH --job-name=umbrella
#SBATCH --partition=gpu
#SBATCH --array=0-25
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=umbrella_%A_%a.out
#SBATCH --error=umbrella_%A_%a.err

module load gromacs/2024.2-cuda

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Window position from array index (0.0 to 5.0 nm, 0.2 nm spacing)
WINDOW=$(echo "scale=1; $SLURM_ARRAY_TASK_ID * 0.2" | bc)

gmx mdrun -deffnm umbrella_${WINDOW} \
    -ntmpi 1 \
    -ntomp $OMP_NUM_THREADS \
    -nb gpu \
    -pme gpu
```

## Continuation Job (Checkpoint Restart)

```bash
#!/bin/bash
#SBATCH --job-name=gromacs_continue
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

module load gromacs/2024.2-cuda

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Continue from checkpoint (appends to existing files)
gmx mdrun -deffnm md \
    -cpi state.cpt \
    -ntmpi 1 \
    -ntomp $OMP_NUM_THREADS \
    -nb gpu \
    -pme gpu \
    -bonded gpu \
    -update gpu
```

## Wall-Time Limited with Auto-Checkpoint

```bash
#!/bin/bash
#SBATCH --job-name=gromacs_walltime
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

module load gromacs/2024.2-cuda

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Automatically stop at 23.5 hours and write checkpoint
# Requires nsteps = -1 in .mdp
gmx mdrun -deffnm md \
    -cpi state.cpt \
    -maxh 23.5 \
    -ntmpi 1 \
    -ntomp $OMP_NUM_THREADS \
    -nb gpu \
    -pme gpu
```

## Complete Workflow Script

```bash
#!/bin/bash
#SBATCH --job-name=gromacs_workflow
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

module load gromacs/2024.2-cuda

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
GMX_OPTS="-ntmpi 1 -ntomp $OMP_NUM_THREADS -nb gpu -pme gpu"

set -e  # Exit on error

echo "=== Energy Minimization ==="
gmx grompp -f em.mdp -c system.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em $GMX_OPTS

# Check EM converged
if ! grep -q "Steepest Descents converged" em.log; then
    echo "ERROR: Energy minimization did not converge"
    exit 1
fi

echo "=== NVT Equilibration ==="
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt $GMX_OPTS

echo "=== NPT Equilibration ==="
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt $GMX_OPTS

echo "=== Production MD ==="
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md $GMX_OPTS -cpo state.cpt

echo "=== Complete ==="
```

## Performance Tuning Tips

### GPU Offloading Flags
```bash
-nb gpu        # Non-bonded on GPU (almost always)
-pme gpu       # PME on GPU (CUDA ≥11)
-bonded gpu    # Bonded forces on GPU (small systems may be faster on CPU)
-update gpu    # Integration on GPU (best for single-GPU)
```

### Optimal Settings by System Size

| System Size | Recommended Setup |
|-------------|-------------------|
| < 50k atoms | 1 GPU, -update gpu |
| 50-200k atoms | 1-2 GPUs |
| 200k-1M atoms | 2-4 GPUs or multi-node |
| > 1M atoms | Multi-node MPI |

### nstlist Tuning
```mdp
; Larger nstlist = better GPU efficiency, but more memory
nstlist = 40    ; Good for GPU (default 10 often too low)
```

### Checking Performance
```bash
# Look at end of log file
tail -50 md.log | grep -A 20 "Performance"

# Key metrics:
# - ns/day (higher is better)
# - Core-hours (lower is better)
# - GPU utilization should be high
```
