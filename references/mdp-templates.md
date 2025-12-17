# MDP Templates Reference

Complete .mdp examples for common simulation stages.

## Energy Minimization

```mdp
; em.mdp - Steepest descent energy minimization
integrator  = steep
emtol       = 1000.0    ; kJ/mol/nm convergence criterion
emstep      = 0.01      ; nm step size
nsteps      = 50000     ; max iterations

; Output control
nstxout     = 500
nstlog      = 500
nstenergy   = 500

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 10
ns_type         = grid
pbc             = xyz
rlist           = 1.2

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.2
pme_order       = 4
fourierspacing  = 0.12

; Van der Waals
vdwtype         = Cut-off
rvdw            = 1.2
DispCorr        = EnerPres

; No velocity generation for EM
gen_vel         = no
```

## Ion Placement (minimal .mdp)

```mdp
; ions.mdp - Used only for grompp before genion
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000

cutoff-scheme   = Verlet
nstlist         = 10
pbc             = xyz
coulombtype     = PME
rcoulomb        = 1.2
rvdw            = 1.2
```

## NVT Equilibration

```mdp
; nvt.mdp - Temperature equilibration
define      = -DPOSRES      ; Position restraints on heavy atoms

integrator  = md
dt          = 0.002         ; 2 fs timestep
nsteps      = 50000         ; 100 ps

; Output
nstxout-compressed  = 5000  ; xtc every 10 ps
nstlog              = 5000
nstenergy           = 5000

; Bond constraints
continuation    = no
constraint_algorithm = LINCS
constraints     = h-bonds

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 20
rlist           = 1.2
pbc             = xyz

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.2

; Van der Waals
vdwtype         = Cut-off
vdw-modifier    = Force-switch
rvdw_switch     = 1.0
rvdw            = 1.2

; Temperature coupling
tcoupl      = V-Rescale
tc_grps     = SYSTEM
tau_t       = 1.0
ref_t       = 300

; No pressure coupling in NVT
pcoupl      = no

; Center of mass motion removal
nstcomm     = 100
comm_mode   = linear
comm_grps   = SYSTEM

; Velocity generation
gen_vel     = yes
gen_temp    = 300
gen_seed    = -1            ; Random seed
```

## NPT Equilibration

```mdp
; npt.mdp - Pressure equilibration
define      = -DPOSRES

integrator  = md
dt          = 0.002
nsteps      = 50000         ; 100 ps

; Output
nstxout-compressed  = 5000
nstlog              = 5000
nstenergy           = 5000

; Continuation from NVT
continuation    = yes
constraint_algorithm = LINCS
constraints     = h-bonds

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 20
rlist           = 1.2
pbc             = xyz

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.2

; Van der Waals
vdwtype         = Cut-off
vdw-modifier    = Force-switch
rvdw_switch     = 1.0
rvdw            = 1.2

; Temperature coupling
tcoupl      = V-Rescale
tc_grps     = SYSTEM
tau_t       = 1.0
ref_t       = 300

; Pressure coupling
pcoupl          = C-rescale
pcoupltype      = isotropic
tau_p           = 5.0
ref_p           = 1.0
compressibility = 4.5e-5

; Center of mass motion removal
nstcomm     = 100
comm_mode   = linear
comm_grps   = SYSTEM

; Velocity generation
gen_vel     = yes
gen_temp    = 300
gen_seed    = -1
```

## Production MD

```mdp
; md.mdp - Production run
integrator  = md
dt          = 0.002
nsteps      = 5000000       ; 10 ns

; Output (balance storage vs resolution)
nstxout-compressed  = 5000  ; 10 ps
nstvout             = 0     ; No velocities saved
nstfout             = 0     ; No forces saved
nstlog              = 10000
nstenergy           = 5000
compressed-x-grps   = System

; Continuation from NPT
continuation    = yes
constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4

; Neighbor searching
cutoff-scheme   = Verlet
nstlist         = 20
pbc             = xyz

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.2
pme_order       = 4
fourierspacing  = 0.12

; Van der Waals
vdwtype         = Cut-off
rvdw            = 1.2
DispCorr        = EnerPres

; Temperature coupling (separate groups often better)
tcoupl      = v-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1   0.1
ref_t       = 300   300

; Pressure coupling
pcoupl          = C-rescale
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = 1.0
compressibility = 4.5e-5

gen_vel     = no
```

## Extended Production (Wall-Time Limited)

```mdp
; md_extended.mdp - For long runs with checkpointing
integrator  = md
dt          = 0.002
nsteps      = -1            ; Run until maxh limit

; ... (same as production above) ...

; Use with: gmx mdrun -maxh 23.5 -cpi state.cpt
```

## High-Frequency Output (Analysis)

```mdp
; md_analysis.mdp - Short run with frequent output for detailed analysis
nstxout-compressed  = 100   ; 0.2 ps resolution
nstenergy           = 100
nstlog              = 1000

; Useful for: MSD fitting, velocity autocorrelation, fast dynamics
```

## Position Restraints Definition

Add to topology or use `#ifdef POSRES`:

```
; In topol.top or posre.itp
#ifdef POSRES
#include "posre.itp"
#endif
```

Generate restraint file:
```bash
gmx genrestr -f em.gro -o posre.itp -fc 1000 1000 1000
```

## Force Field Specific Notes

### CHARMM36
- Use `charmm36-jul2022.ff` or newer
- Water: `tip3p` (modified TIP3P)
- Cutoffs: 1.2 nm (force-switch for vdW recommended)
- Add: `vdw-modifier = Force-switch`, `rvdw-switch = 1.0`

### AMBER (ff19SB)
- Water: `tip3p` or `opc`
- Cutoffs: 0.9-1.0 nm typical
- PME for electrostatics

### OPLS-AA
- Water: `tip4p` or `tip3p`
- Cutoffs: 1.0-1.2 nm
