# Free Energy Methods Reference

Guide for umbrella sampling, steered MD, PMF calculations, and alchemical transformations.

## Umbrella Sampling Workflow

### 1. Generate Initial Configurations

Use steered MD (pulling) to generate windows along reaction coordinate:

```bash
# Pull simulation
gmx grompp -f pull.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o pull.tpr
gmx mdrun -deffnm pull -pf pullf.xvg -px pullx.xvg
```

### 2. Pull MDP Parameters

```mdp
; pull.mdp - Steered MD for umbrella sampling preparation
integrator  = md
dt          = 0.002
nsteps      = 500000        ; 1 ns pulling

; Output
nstxout             = 5000  ; Save conformations for window extraction
nstxout-compressed  = 500
nstlog              = 500
nstenergy           = 500

; ... standard MD parameters ...

; Pull code
pull                    = yes
pull-ncoords            = 1
pull-ngroups            = 2
pull-group1-name        = Ligand
pull-group2-name        = Protein

pull-coord1-type        = umbrella
pull-coord1-geometry    = distance
pull-coord1-dim         = Y Y Y
pull-coord1-groups      = 1 2
pull-coord1-start       = yes
pull-coord1-rate        = 0.01      ; nm/ps (10 nm/ns)
pull-coord1-k           = 1000      ; kJ/mol/nm² force constant
```

### 3. Extract Window Configurations

```bash
# Extract frames at regular intervals along pull coordinate
for i in $(seq 0 0.2 5.0); do
    gmx trjconv -f pull.xtc -s pull.tpr -o conf_${i}.gro \
        -dump $(echo "$i * 100" | bc)  # Time depends on pull rate
done
```

Or use `gmx select` with distance criteria.

### 4. Umbrella Window Simulations

```mdp
; umbrella.mdp - Window simulation
integrator  = md
dt          = 0.002
nsteps      = 5000000       ; 10 ns per window

; Standard parameters...

; Umbrella restraint
pull                    = yes
pull-ncoords            = 1
pull-ngroups            = 2
pull-group1-name        = Ligand
pull-group2-name        = Protein

pull-coord1-type        = umbrella
pull-coord1-geometry    = distance
pull-coord1-dim         = Y Y Y
pull-coord1-groups      = 1 2
pull-coord1-init        = WINDOW_POSITION  ; Replace per window
pull-coord1-k           = 1000             ; Same as pulling
pull-coord1-rate        = 0                ; Static umbrella

; Output pull data
pull-nstxout            = 100
pull-nstfout            = 100
```

Run each window:
```bash
for i in $(seq 0 0.2 5.0); do
    sed "s/WINDOW_POSITION/${i}/" umbrella_template.mdp > umbrella_${i}.mdp
    gmx grompp -f umbrella_${i}.mdp -c conf_${i}.gro -p topol.top -n index.ndx -o umbrella_${i}.tpr
    gmx mdrun -deffnm umbrella_${i} -pf pullf_${i}.xvg -px pullx_${i}.xvg
done
```

### 5. WHAM Analysis

```bash
# Create file lists
ls pullf_*.xvg > files_pullf.dat
ls umbrella_*.tpr > files_tpr.dat

# Run WHAM
gmx wham -it files_tpr.dat -if files_pullf.dat \
    -o pmf.xvg -hist hist.xvg -bsres bsResult.xvg \
    -nBootstrap 100 -bs-method hist

# Options:
# -temp 300          Temperature
# -bins 100          Number of bins
# -min 0 -max 5      Range (nm)
# -tol 1e-6          Convergence tolerance
```

### 6. Convergence Analysis

```bash
# Time-dependent PMF (check convergence)
for t in 2 4 6 8 10; do
    gmx wham -it files_tpr.dat -if files_pullf.dat \
        -o pmf_${t}ns.xvg -b $((t-2))000 -e ${t}000
done
```

## Steered MD Variants

### Constant Velocity (SMD)
```mdp
pull-coord1-type    = umbrella
pull-coord1-rate    = 0.001     ; nm/ps
pull-coord1-k       = 1000      ; kJ/mol/nm²
```

### Constant Force
```mdp
pull-coord1-type    = constant-force
pull-coord1-k       = 500       ; kJ/mol/nm (force, not spring constant)
```

## Free Energy Perturbation (FEP)

### Lambda Schedule

```mdp
; fep.mdp - Alchemical transformation
free-energy         = yes
init-lambda-state   = 0         ; Change per window (0 to n-1)

; Soft-core parameters (avoid singularities)
sc-alpha            = 0.5
sc-power            = 1
sc-sigma            = 0.3

; Lambda vectors (20 windows example)
; vdw-lambdas: Turn off vdW first
; coul-lambdas: Then electrostatics
vdw-lambdas         = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
coul-lambdas        = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
bonded-lambdas      = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0

; Output dH/dlambda
nstdhdl             = 100
calc-lambda-neighbors = -1      ; All lambdas for MBAR
```

### Analysis with BAR/MBAR

```bash
# Using gmx bar
gmx bar -f dhdl_*.xvg -o bar.xvg -oi barint.xvg

# Or use alchemlyb (Python) for MBAR
pip install alchemlyb
```

```python
# alchemlyb analysis
from alchemlyb.parsing.gmx import extract_dHdl, extract_u_nk
from alchemlyb.estimators import MBAR, BAR

u_nk = pd.concat([extract_u_nk(f, T=300) for f in sorted(glob('dhdl_*.xvg'))])
mbar = MBAR().fit(u_nk)
print(f"ΔG = {mbar.delta_f_.iloc[0, -1]:.2f} kJ/mol")
```

## Common Issues

### Umbrella Sampling
- **Poor overlap**: Reduce window spacing or increase force constant
- **Insufficient sampling**: Extend simulations, check autocorrelation
- **Asymmetric PMF**: Check for hysteresis, may need more equilibration

### FEP
- **Soft-core issues**: Adjust sc-alpha and sc-sigma
- **Poor convergence**: Add intermediate lambda states
- **End-state instabilities**: Use longer equilibration at λ=0 and λ=1

## Analysis Scripts

### Check window histogram overlap
```python
import numpy as np
import matplotlib.pyplot as plt

# Load WHAM histograms
data = np.loadtxt('hist.xvg', comments=['#', '@'])
plt.figure(figsize=(10, 6))
for i in range(1, data.shape[1]):
    plt.plot(data[:, 0], data[:, i])
plt.xlabel('Distance (nm)')
plt.ylabel('Count')
plt.title('Umbrella Sampling Histogram Overlap')
plt.savefig('histogram_overlap.png', dpi=150)
```

### PMF convergence plot
```python
import numpy as np
import matplotlib.pyplot as plt

times = [2, 4, 6, 8, 10]
plt.figure(figsize=(10, 6))
for t in times:
    data = np.loadtxt(f'pmf_{t}ns.xvg', comments=['#', '@'])
    plt.plot(data[:, 0], data[:, 1], label=f'{t} ns')
plt.xlabel('Distance (nm)')
plt.ylabel('PMF (kJ/mol)')
plt.legend()
plt.title('PMF Convergence')
plt.savefig('pmf_convergence.png', dpi=150)
```
