# Swarm HPC Workflow

Guide for preparing GROMACS simulations on workstation and submitting to swarm cluster.

## Environment Configuration

### Workstation (Local)
- **GROMACS command**: `gmx_mpi`
- **Purpose**: File preparation, topology generation, preprocessing
- **No SLURM**: Run commands directly

### Swarm Cluster (swarm.che.wisc.edu)
- **GROMACS command**: `gmx`
- **Purpose**: Production MD runs
- **SLURM**: Use job scheduler
- **Access**: `ssh aaltamimi2@swarm.che.wisc.edu`

## Standard Workflow

### Step 1: Prepare Files on Workstation

```bash
# 1. Generate topology (workstation uses gmx_mpi)
gmx_mpi pdb2gmx -f protein.pdb -o protein.gro -p topol.top -ff amber99sb-ildn

# 2. Define box and solvate
gmx_mpi editconf -f protein.gro -o boxed.gro -c -d 1.2 -bt dodecahedron
gmx_mpi solvate -cp boxed.gro -o solvated.gro -p topol.top

# 3. Add ions
gmx_mpi grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
echo "SOL" | gmx_mpi genion -s ions.tpr -o system.gro -p topol.top -neutral

# 4. Energy minimization (preprocessing only - check system)
gmx_mpi grompp -f em.mdp -c system.gro -p topol.top -o em.tpr
```

### Step 2: Create SLURM Job Script

**Template**: `submit_swarm.sh`

```bash
#!/bin/bash
#SBATCH -p compute
#SBATCH -t 108:00:00
#SBATCH -J {JOB_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mail-user=aaltamimi2@wisc.edu
#SBATCH --mail-type=end

# Load GROMACS module on swarm
module load gromacs

# Energy Minimization
gmx mdrun -v -deffnm em

# NVT Equilibration
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

# NPT Equilibration
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt

# Production MD
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md
```

**Variable Substitution**:
- `{JOB_NAME}`: Replace with your job name (e.g., `protein_md`, `ligand_binding`)

### Step 3: Transfer Files to Swarm

```bash
# Create project directory on swarm
ssh aaltamimi2@swarm.che.wisc.edu "mkdir -p ~/projects/my_simulation"

# Transfer input files
scp system.gro topol.top *.mdp *.itp submit_swarm.sh \
    aaltamimi2@swarm.che.wisc.edu:~/projects/my_simulation/

# If using ligands, transfer ACPYPE outputs
scp -r LIG.acpype aaltamimi2@swarm.che.wisc.edu:~/projects/my_simulation/
```

### Step 4: Submit Job on Swarm

```bash
# SSH to swarm
ssh aaltamimi2@swarm.che.wisc.edu

# Navigate to project
cd ~/projects/my_simulation

# Submit job
sbatch submit_swarm.sh

# Check job status
squeue -u aaltamimi2

# Monitor output
tail -f slurm-*.out
```

### Step 5: Retrieve Results

```bash
# From workstation, download results
scp -r aaltamimi2@swarm.che.wisc.edu:~/projects/my_simulation/*.xtc .
scp -r aaltamimi2@swarm.che.wisc.edu:~/projects/my_simulation/*.edr .
scp -r aaltamimi2@swarm.che.wisc.edu:~/projects/my_simulation/*.log .
```

## SLURM Job Templates

### Template 1: Energy Minimization + Equilibration + Production

```bash
#!/bin/bash
#SBATCH -p compute
#SBATCH -t 108:00:00
#SBATCH -J protein_md
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mail-user=aaltamimi2@wisc.edu
#SBATCH --mail-type=end

module load gromacs

# EM
gmx mdrun -v -deffnm em

# NVT
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

# NPT
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt

# Production
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md
```

### Template 2: Production Only (Equilibration Done)

```bash
#!/bin/bash
#SBATCH -p compute
#SBATCH -t 108:00:00
#SBATCH -J md_production
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mail-user=aaltamimi2@wisc.edu
#SBATCH --mail-type=end

module load gromacs

# Production MD
gmx mdrun -deffnm md
```

### Template 3: Continuation Run

```bash
#!/bin/bash
#SBATCH -p compute
#SBATCH -t 108:00:00
#SBATCH -J md_continue
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mail-user=aaltamimi2@wisc.edu
#SBATCH --mail-type=end

module load gromacs

# Continue from checkpoint
gmx mdrun -deffnm md -cpi state.cpt -append
```

### Template 4: Array Job (Multiple Simulations)

```bash
#!/bin/bash
#SBATCH -p compute
#SBATCH -t 108:00:00
#SBATCH -J replica_array
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --array=1-10
#SBATCH --mail-user=aaltamimi2@wisc.edu
#SBATCH --mail-type=end

module load gromacs

# Create replica directory
REPLICA_DIR="replica_${SLURM_ARRAY_TASK_ID}"
mkdir -p $REPLICA_DIR
cd $REPLICA_DIR

# Copy input files
cp ../em.tpr ../nvt.mdp ../npt.mdp ../md.mdp ../topol.top .

# Run with different seed
gmx mdrun -deffnm em
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md
```

## Helper Commands

### SSH Key Setup (One-Time)

```bash
# Generate SSH key on workstation
ssh-keygen -t rsa -b 4096

# Copy to swarm for passwordless login
ssh-copy-id aaltamimi2@swarm.che.wisc.edu
```

### Quick Transfer Script

Create `transfer_to_swarm.sh` on workstation:

```bash
#!/bin/bash
# Quick transfer to swarm

PROJECT_NAME=$1
if [[ -z "$PROJECT_NAME" ]]; then
    echo "Usage: $0 <project_name>"
    exit 1
fi

SWARM_DIR="~/projects/${PROJECT_NAME}"

# Create directory on swarm
ssh aaltamimi2@swarm.che.wisc.edu "mkdir -p $SWARM_DIR"

# Transfer all necessary files
scp *.gro *.top *.mdp *.itp submit_swarm.sh \
    aaltamimi2@swarm.che.wisc.edu:$SWARM_DIR/

# Transfer ligand directories if they exist
for dir in *.acpype; do
    if [[ -d "$dir" ]]; then
        scp -r "$dir" aaltamimi2@swarm.che.wisc.edu:$SWARM_DIR/
    fi
done

echo "Files transferred to swarm:$SWARM_DIR"
```

Usage:
```bash
chmod +x transfer_to_swarm.sh
./transfer_to_swarm.sh protein_ligand_md
```

## Job Management

### Check Job Status

```bash
# View your jobs
squeue -u aaltamimi2

# Detailed job info
scontrol show job <JOB_ID>

# View job history
sacct -u aaltamimi2

# Cancel job
scancel <JOB_ID>
```

### Monitor Running Job

```bash
# Watch output file
tail -f slurm-<JOB_ID>.out

# Check log file
tail -f md.log

# Monitor in real-time
watch -n 10 'tail -20 md.log'
```

## Important Notes

### Command Differences

| Location | Command | Example |
|----------|---------|---------|
| **Workstation** | `gmx_mpi` | `gmx_mpi grompp -f md.mdp ...` |
| **Swarm** | `gmx` | `gmx mdrun -deffnm md` |

### Time Limits

- Maximum time: 108 hours (4.5 days)
- For longer runs: Use checkpoints and continuation jobs
- Plan for automatic restarts

### File Organization

```
project_directory/
├── input/
│   ├── protein.pdb
│   ├── ligand.mol2
│   └── *.mdp
├── topology/
│   ├── topol.top
│   ├── *.itp
│   └── LIG.acpype/
├── submit_swarm.sh
└── analysis/
    └── (created after job completes)
```

## Troubleshooting

### Cannot Connect to Swarm

```bash
# Test connection
ssh -v aaltamimi2@swarm.che.wisc.edu

# Check VPN if off-campus
# University of Wisconsin VPN required
```

### Job Fails Immediately

```bash
# Check SLURM output
cat slurm-*.out

# Verify module loads
module avail gromacs
module load gromacs
which gmx
```

### Files Not Found on Swarm

```bash
# Verify transfer
ssh aaltamimi2@swarm.che.wisc.edu "ls -la ~/projects/my_simulation"

# Check paths in job script
# All paths should be relative or explicit
```

## Best Practices

1. **Test Locally First**: Run short test on workstation with `gmx_mpi` before submitting to swarm
2. **Use Descriptive Job Names**: `{JOB_NAME}` should indicate project, system, and purpose
3. **Checkpoint Frequently**: Set `nstxout-compressed` appropriately in .mdp
4. **Email Notifications**: Enabled by `#SBATCH --mail-type=end`
5. **Backup**: Keep local copies of all input files before transferring
