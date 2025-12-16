# GROMACS Session Hook Example

Create `~/.claude/hooks/session-start.sh` to automatically set up your GROMACS environment:

```bash
#!/bin/bash
# ~/.claude/hooks/session-start.sh

# Load GROMACS module if on HPC
if command -v module &> /dev/null; then
    module load gromacs/2021.7 2>/dev/null || true
fi

# Set up GROMACS environment
export GMXLIB=/home/aaltamimi2/gromacs-plumed/share/gromacs/top
export GMX=/home/aaltamimi2/gromacs-plumed/bin/gmx_mpi

# Create alias for MPI wrapper
alias gmx='mpirun -np 1 $GMX'

# Print GROMACS version for confirmation
echo "GROMACS environment loaded:"
$GMX --version 2>&1 | head -1
echo "Use 'gmx' command (auto-wraps with mpirun)"
```

Make it executable:
```bash
chmod +x ~/.claude/hooks/session-start.sh
```

Now every Claude Code session will have GROMACS ready!
