# Expanding Your Claude Code Skills Usage

A guide to getting the most out of Claude Code skills for GROMACS and beyond.

## 🚀 Quick Wins (Start Here)

### 1. Use the Skill in Conversations

Simply ask Claude Code for GROMACS help:

```
"Help me set up a membrane protein simulation in POPC lipid bilayer"
"Create production .mdp for 500ns NPT at 310K"
"My simulation has LINCS warnings - help me debug"
"Analyze my trajectory for protein-ligand binding stability"
```

Claude will automatically use the GROMACS skill knowledge.

### 2. Create Project-Specific Commands

```bash
# Create commands directory
mkdir -p ~/.claude/commands

# Example: Quick validation command
cat > ~/.claude/commands/gmx-validate.md << 'EOF'
---
description: "Validate GROMACS installation"
---
Run the GROMACS validation test and verify the installation is working.
Use the validation script and check for temperature ~300K.
EOF
```

Use with: `/gmx-validate` in Claude Code

### 3. Set Up Environment Alias

Add to your `~/.bashrc`:

```bash
# GROMACS shortcuts
export GMX=/home/aaltamimi2/gromacs-plumed/bin/gmx_mpi
alias gmx='mpirun -np 1 $GMX'
alias gmx-help='claude "GROMACS help:"'
```

## 📚 Intermediate: Customize for Your Research

### Create Research-Specific Skills

Structure:
```
~/.claude/skills/
├── gromacs-skill/          # Base GROMACS skill (from repo)
├── membrane-proteins/      # Your custom skill
│   ├── SKILL.md
│   └── references/
│       ├── popc-bilayer.md
│       ├── gpcr-protocols.md
│       └── apo-parameters.md
└── free-energy/            # Another custom skill
    ├── SKILL.md
    └── references/
        ├── umbrella-sampling.md
        └── alchemical-fep.md
```

Example custom skill (`~/.claude/skills/membrane-proteins/SKILL.md`):

```markdown
---
name: membrane-proteins
description: "Specialized protocols for membrane protein MD simulations with GROMACS"
---

# Membrane Protein Simulations

## System Building with CHARMM-GUI

1. Upload PDB to CHARMM-GUI Membrane Builder
2. Select membrane type (POPC, POPE, etc.)
3. Download GROMACS files
4. Continue with equilibration

## Equilibration Protocol

6-step equilibration (CHARMM-GUI default):
```bash
for i in {1..6}; do
  mpirun -np 1 gmx_mpi grompp -f step6.${i}_equilibration.mdp -c step6.${i}.gro -p topol.top -o step6.${i}.tpr
  mpirun -np 8 gmx_mpi mdrun -deffnm step6.${i}
done
```

## Production MD

Use NPT with semi-isotropic pressure coupling for membranes:
```mdp
pcoupl          = C-rescale
pcoupltype      = semiisotropic
tau_p           = 5.0 5.0
ref_p           = 1.0 1.0
compressibility = 4.5e-5 4.5e-5
```

## Common Issues

- **Membrane thinning**: Increase tau_p to 5.0-10.0
- **Protein tilting**: Check initial orientation, use position restraints longer
- **Lipid scrambling**: Reduce time step to 1 fs initially
```

### Create Analysis Templates Library

```bash
mkdir -p ~/gromacs-templates/analysis

# Save reusable analysis scripts
cp ~/gromacs-claudecode-skill/examples/analysis-workflows.md ~/gromacs-templates/
```

Ask Claude: "Use the standard protein analysis workflow from my templates"

## 🔥 Advanced: Full Integration

### 1. Automated Workflows with Scripts

Create `~/bin/gromacs-workflow.sh`:

```bash
#!/bin/bash
# Automated GROMACS workflow with Claude Code assistance

SYSTEM=$1
if [[ -z "$SYSTEM" ]]; then
  echo "Usage: $0 <system_name>"
  exit 1
fi

# Ask Claude to generate the workflow
claude << EOF
Create a complete GROMACS workflow for $SYSTEM including:
1. System preparation (pdb2gmx, solvate, ions)
2. Energy minimization
3. NVT equilibration
4. NPT equilibration
5. Production MD (100 ns)
6. Analysis pipeline

Generate all .mdp files and shell scripts.
EOF
```

### 2. Git-Tracked Simulation Projects

```bash
# Template repository structure
my-simulation-project/
├── .claude/
│   └── commands/
│       ├── setup-system.md
│       ├── run-equilibration.md
│       └── analyze-trajectory.md
├── inputs/
│   └── protein.pdb
├── parameters/
│   ├── em.mdp
│   ├── nvt.mdp
│   ├── npt.mdp
│   └── production.mdp
├── scripts/
│   ├── 01-setup.sh
│   ├── 02-minimize.sh
│   ├── 03-equilibrate.sh
│   ├── 04-production.sh
│   └── 05-analyze.sh
└── README.md
```

Ask Claude to help maintain this structure for each project.

### 3. MCP Servers for GROMACS

Create a Model Context Protocol server for GROMACS:

```python
# ~/.claude/mcp-servers/gromacs/server.py
# Provides GROMACS commands as MCP tools

from mcp.server import Server
import subprocess

server = Server("gromacs-mcp")

@server.tool("run_gmx_command")
def run_gmx(command: str, args: list[str]) -> str:
    """Execute GROMACS command via MPI wrapper"""
    cmd = ["mpirun", "-np", "1", "gmx_mpi", command] + args
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout

@server.tool("check_simulation_status")
def check_status(log_file: str) -> dict:
    """Parse GROMACS log file for status"""
    # Parse em.log, nvt.log, etc.
    # Return: {converged: bool, energy: float, time: float}
    pass
```

### 4. Custom Slash Commands for Common Tasks

```bash
# ~/.claude/commands/setup-membrane.md
---
description: "Set up membrane protein system with CHARMM-GUI files"
---
Guide me through setting up a membrane protein simulation using CHARMM-GUI output files.
Include proper equilibration protocol and production MD setup.
```

```bash
# ~/.claude/commands/analyze-binding.md
---
description: "Analyze protein-ligand binding from trajectory"
---
Run comprehensive protein-ligand binding analysis including:
- Ligand RMSD
- Protein-ligand distances
- H-bond analysis
- Binding pocket RMSF
- Generate publication-quality plots
```

## 📊 Integration with HPC

### SLURM Job Templates

Ask Claude to generate job scripts:

```
"Create a SLURM job script for 1TB node with 4 A100 GPUs running
500ns GROMACS production MD. Use optimal GPU offloading settings."
```

### Module Integration

```bash
# ~/.claude/hooks/session-start.sh
#!/bin/bash

# Auto-load modules on HPC
if [[ -f /etc/profile.d/modules.sh ]]; then
  source /etc/profile.d/modules.sh
  module load gromacs/2023.1-gpu
  module load cuda/11.8
  module load openmpi/4.1.4
fi

# Set up environment
export GMXLIB=$EBROOTGROMACS/share/gromacs/top
alias gmx='mpirun -np 1 gmx_mpi'

echo "✓ GROMACS environment ready (HPC mode)"
```

## 🎯 Real-World Usage Examples

### Example 1: New Simulation Setup

```bash
# In your project directory
claude "Set up a new GROMACS simulation for protein.pdb:
1. Use AMBER99SB-ILDN force field
2. TIP3P water in dodecahedron box
3. Neutralize and add 150mM NaCl
4. Generate all necessary .mdp files
5. Create run scripts for local testing and HPC submission"
```

### Example 2: Troubleshooting

```bash
# When you hit an error
claude "My simulation crashed with LINCS warning. Here's my .mdp:
[paste em.mdp]

And the error:
[paste error message]

Help me fix this and explain what went wrong."
```

### Example 3: Analysis Pipeline

```bash
claude "Analyze md.xtc for:
1. Protein RMSD and RMSF
2. Secondary structure evolution
3. Key H-bonds stability
4. Generate plots with matplotlib
Create a complete analysis script."
```

## 🔧 Maintenance & Updates

### Keep Skills Updated

```bash
cd ~/gromacs-claudecode-skill
git pull origin main

# Skills in .claude/skills/ update automatically via symlinks
```

### Share Skills with Team

```bash
# Fork the repository
# Add your custom references
git add references/my-custom-protocol.md
git commit -m "Add custom protocol for X"
git push

# Team members can pull your updates
```

## 📖 Resources

- GROMACS Manual: http://manual.gromacs.org/
- Claude Code Docs: Check skill README
- Your custom skills: `~/.claude/skills/`
- Templates: `~/gromacs-templates/`

## Next Steps

1. ✅ Install skill (done - you have this working!)
2. ✅ Test with validation script (done!)
3. ⬜ Create first custom command (`/gmx-validate`)
4. ⬜ Set up analysis templates
5. ⬜ Create research-specific skill extension
6. ⬜ Integrate with HPC workflow
7. ⬜ Share with lab members

Start with steps 3-4 this week, then expand gradually!
