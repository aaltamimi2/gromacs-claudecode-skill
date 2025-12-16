#!/bin/bash
# Generate SLURM job script for swarm cluster
# Uses EXACTLY the format required for swarm.che.wisc.edu

usage() {
    cat << EOF
Usage: $0 -n JOB_NAME [-t TYPE] [-o OUTPUT]

Generate SLURM job script for swarm cluster (swarm.che.wisc.edu)

Options:
    -n JOB_NAME     Job name (required)
    -t TYPE         Job type: full, production, continue, array (default: full)
    -o OUTPUT       Output filename (default: submit_swarm.sh)
    -h              Show this help message

Job Types:
    full        - Complete workflow: EM + NVT + NPT + Production
    production  - Production MD only (equilibration already done)
    continue    - Continue from checkpoint
    array       - Array job for multiple replicas

Examples:
    $0 -n protein_md -t full
    $0 -n md_prod -t production -o submit_production.sh
    $0 -n replica_1 -t continue
EOF
    exit 1
}

# Default values
JOB_TYPE="full"
OUTPUT="submit_swarm.sh"

# Parse arguments
while getopts "n:t:o:h" opt; do
    case $opt in
        n) JOB_NAME="$OPTARG" ;;
        t) JOB_TYPE="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check required arguments
if [[ -z "$JOB_NAME" ]]; then
    echo "Error: Job name (-n) is required"
    usage
fi

# Validate job type
if [[ ! "$JOB_TYPE" =~ ^(full|production|continue|array)$ ]]; then
    echo "Error: Invalid job type. Must be: full, production, continue, or array"
    usage
fi

# Generate job script based on type
cat > "$OUTPUT" << 'HEADER'
#!/bin/bash
#SBATCH -p compute
#SBATCH -t 108:00:00
HEADER

echo "#SBATCH -J $JOB_NAME" >> "$OUTPUT"

cat >> "$OUTPUT" << 'FOOTER'
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mail-user=aaltamimi2@wisc.edu
#SBATCH --mail-type=end

module load gromacs

FOOTER

# Add job-specific commands
case $JOB_TYPE in
    full)
        cat >> "$OUTPUT" << 'EOF'
# Energy Minimization
gmx mdrun -nt 28 -v -deffnm em

# NVT Equilibration
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt -nb gpu -bonded cpu -update gpu -ntomp 16 -pin on

# NPT Equilibration
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt -nb gpu -bonded cpu -update gpu -ntomp 16 -pin on

# Production MD
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -v -deffnm md -nb gpu -bonded cpu -update gpu -ntomp 16 -pin on
EOF
        ;;

    production)
        cat >> "$OUTPUT" << 'EOF'
# Production MD (equilibration already done)
gmx mdrun -v -deffnm md -nb gpu -bonded cpu -update gpu -ntomp 16 -pin on
EOF
        ;;

    continue)
        cat >> "$OUTPUT" << 'EOF'
# Continue simulation from checkpoint
gmx mdrun -v -deffnm md -cpi state.cpt -append -nb gpu -bonded cpu -update gpu -ntomp 16 -pin on
EOF
        ;;

    array)
        # Update SBATCH header for array job
        sed -i "s/#SBATCH --mail-type=end/#SBATCH --array=1-10\n#SBATCH --mail-type=end/" "$OUTPUT"

        cat >> "$OUTPUT" << 'EOF'
# Create replica directory
REPLICA_DIR="replica_${SLURM_ARRAY_TASK_ID}"
mkdir -p $REPLICA_DIR
cd $REPLICA_DIR

# Copy input files
cp ../em.tpr ../nvt.mdp ../npt.mdp ../md.mdp ../topol.top .

# Run simulation
gmx mdrun -nt 28 -v -deffnm em
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt -nb gpu -bonded cpu -update gpu -ntomp 16 -pin on

gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt -nb gpu -bonded cpu -update gpu -ntomp 16 -pin on

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -v -deffnm md -nb gpu -bonded cpu -update gpu -ntomp 16 -pin on
EOF
        ;;
esac

chmod +x "$OUTPUT"

echo "✓ Generated SLURM job script: $OUTPUT"
echo "  Job name: $JOB_NAME"
echo "  Job type: $JOB_TYPE"
echo ""
echo "Next steps:"
echo "  1. Review the script: cat $OUTPUT"
echo "  2. Transfer to swarm: scp $OUTPUT swarm:~/projects/"
echo "  3. Submit on swarm: ssh swarm 'cd ~/projects && sbatch $OUTPUT'"
