# Polymer-Solvent MD Workflow

End-to-end workflow for polymer oligomer MD simulations in organic solvents, designed for COSMO-RS conformational sampling.

## Quick Start

```bash
# Example: PVDF 4mer in NMP
python polymer_md_workflow.py \
  --polymer PVDF \
  --solvent NMP \
  --work-dir ./pvdf_nmp
```

## Workflow Overview

This workflow automates:

1. **Polymer oligomer generation** (with backbone atom validation)
2. **Solvent SMILES validation** (with user confirmation)
3. **Parameterization** (using `solvent_to_gmx.py`)
4. **Solvated box creation** (15,000 solvent molecules in 10×10×10 nm box)
5. **Energy minimization test**
6. **Submission script generation** (swarm cluster and workstation)

### Default Parameters

- **Box size**: 10×10×10 nm
- **Solvent molecules**: 15,000
- **Backbone atoms**: 10-12 (ideal), max 20
- **MD Protocol**:
  - Energy Minimization (EM)
  - 10 ns NVT equilibration
  - 5 ns NPT equilibration

## Usage

### Basic Usage

```bash
python polymer_md_workflow.py --polymer <NAME> --solvent <NAME> --work-dir <DIR>
```

### Options

```
--polymer POLYMER          Polymer name or SMILES (required)
--solvent SOLVENT          Solvent name or SMILES (required)
--work-dir WORK_DIR        Working directory (required)
--oligomer-size N          Number of monomer units (auto if not specified)
--job-name NAME            SLURM job name (auto-generated)
--email EMAIL              Email for notifications (default: aaltamimi2@wisc.edu)
--box-size SIZE            Box size in nm (default: 10.0)
--n-solvent N              Number of solvent molecules (default: 15000)
--target TARGET            Platform: swarm/workstation/both (default: both)
--no-interactive           Run without user prompts
```

### Examples

#### Example 1: PVDF in NMP (Auto Oligomer Size)
```bash
python polymer_md_workflow.py \
  --polymer PVDF \
  --solvent NMP \
  --work-dir ./pvdf_nmp
```

This will:
- Auto-calculate optimal oligomer size to hit ~11 backbone atoms
- Validate NMP SMILES
- Create solvated box with 15,000 NMP molecules
- Generate submission scripts

#### Example 2: Explicit Oligomer Size
```bash
python polymer_md_workflow.py \
  --polymer PVDF \
  --solvent NMP \
  --oligomer-size 4 \
  --work-dir ./pvdf4_nmp
```

Forces a 4-mer (fails if backbone atoms > 20).

#### Example 3: Custom Polymer SMILES
```bash
python polymer_md_workflow.py \
  --polymer "C(C(F)(F))([H])[H]" \
  --solvent "CN1CCCC1=O" \
  --work-dir ./custom_polymer
```

#### Example 4: Non-Interactive Mode
```bash
python polymer_md_workflow.py \
  --polymer PVDF \
  --solvent NMP \
  --work-dir ./test \
  --no-interactive
```

Runs without user confirmation prompts (useful for automation).

## Supported Polymers

| Name | Full Name | SMILES | Good Solvents |
|------|-----------|--------|---------------|
| PVDF | Polyvinylidene fluoride | `C(C(F)(F))([H])[H]` | NMP, DMF, DMSO, acetone |
| PEO  | Polyethylene oxide | `CCO` | water, chloroform, acetonitrile |
| PMMA | Polymethyl methacrylate | `CC(C)(C(=O)OC)` | chloroform, THF, toluene |
| PS   | Polystyrene | `C(C(c1ccccc1))([H])[H]` | toluene, chloroform, THF |

You can also provide custom SMILES strings.

## Supported Solvents

| Name | Full Name | SMILES |
|------|-----------|--------|
| NMP  | N-Methyl-2-pyrrolidone | `CN1CCCC1=O` |
| DMF  | Dimethylformamide | `CN(C)C=O` |
| DMSO | Dimethyl sulfoxide | `CS(=O)C` |
| water | Water | `O` |
| chloroform | Chloroform | `ClC(Cl)Cl` |
| acetone | Acetone | `CC(=O)C` |
| acetonitrile | Acetonitrile | `CC#N` |
| THF  | Tetrahydrofuran | `C1CCOC1` |
| toluene | Toluene | `Cc1ccccc1` |

Custom SMILES are also supported (will be validated).

## Output Files

After successful completion, the working directory contains:

```
work_dir/
├── polymer.acpype/              # Polymer ACPYPE parameterization
├── solvent.acpype/              # Solvent ACPYPE parameterization
├── system_solvated.gro          # Solvated structure
├── topol.top                    # Combined topology
├── em.mdp                       # Energy minimization parameters
├── nvt.mdp                      # NVT parameters (10 ns)
├── npt.mdp                      # NPT parameters (5 ns)
├── submit_swarm.sh              # Swarm cluster submission script
├── run_workstation.sh           # Workstation execution script
└── prepare_simulations.sh       # Prepare .tpr files
```

## Running Simulations

### On Workstation

```bash
cd work_dir
./prepare_simulations.sh      # Prepare .tpr files
./run_workstation.sh           # Run locally (uses gmx_mpi)
```

### On Swarm Cluster

```bash
# 1. Prepare files on workstation
cd work_dir
./prepare_simulations.sh

# 2. Transfer to swarm
scp *.tpr *.mdp *.top *.gro *.itp *.acpype submit_swarm.sh swarm:~/project_dir/

# 3. Submit job
ssh swarm 'cd ~/project_dir && sbatch submit_swarm.sh'

# 4. Monitor
ssh swarm 'squeue -u aaltamimi2'
ssh swarm 'tail -f ~/project_dir/slurm-*.out'
```

## Module Structure

The workflow is organized into modular components:

```
polymer-cosmo/
├── polymer_md_workflow.py          # Main orchestrator
├── modules/
│   ├── polymer_builder.py          # Oligomer generation
│   ├── box_builder.py              # Solvated box creation
│   └── submission_generator.py     # Script generation
├── templates/
│   ├── em.mdp                      # EM parameters
│   ├── nvt.mdp                     # NVT parameters
│   └── npt.mdp                     # NPT parameters
└── README.md
```

### Running Individual Modules

Each module can be run independently for testing/debugging:

#### Test Polymer Builder
```bash
python modules/polymer_builder.py PVDF --recommend-solvents
```

#### Test Box Builder
```bash
python modules/box_builder.py ./work_dir \
  --polymer-smiles "C(C(F)(F))([H])[H]" \
  --solvent-smiles "CN1CCCC1=O" \
  --n-solvent 15000
```

#### Test Submission Generator
```bash
python modules/submission_generator.py ./work_dir \
  --job-name test_job \
  --target both
```

## Integration with GROMACS Skill

This workflow integrates with existing scripts:

- **`scripts/solvent_to_gmx.py`**: Parameterizes polymer and solvent
- **`references/mdp-templates.md`**: MDP parameter documentation
- **`references/swarm-workflow.md`**: Swarm cluster guidelines

## Backbone Atom Constraints

The workflow enforces backbone atom limits:

- **Ideal range**: 10-12 backbone atoms
- **Maximum**: 20 backbone atoms

This ensures:
- Reasonable oligomer size for MD
- Manageable conformational sampling
- Appropriate for DFT calculations

If you specify an oligomer size that exceeds the backbone limit, the workflow will fail with an error.

## Troubleshooting

### "Invalid polymer SMILES"
- Check polymer name spelling (case-insensitive)
- If using custom SMILES, verify syntax with RDKit

### "Solvent insertion failed"
- Box may be too small for 15,000 molecules
- Try reducing `--n-solvent` or increasing `--box-size`

### "Energy minimization test failed"
- Check `em_test.log` in working directory
- May indicate topology issues
- Workflow allows continuing despite EM failure (with confirmation)

### "ACPYPE parameterization failed"
- Ensure RDKit and ACPYPE are installed
- Check SMILES syntax
- Some complex molecules may not be supported by ACPYPE

## Next Steps: Conformational Sampling

After MD simulation completes, the next workflow steps will be:

1. **Conformational sampling** (Rg-SASA grid)
2. **DFT preparation** (ORCA COSMO-RS inputs)

These modules are under development.

## Requirements

- Python 3.7+
- RDKit (`conda install -c conda-forge rdkit`)
- GROMACS (gmx_mpi for workstation, gmx for swarm)
- ACPYPE (installed via `solvent_to_gmx.py`)

## Support

For issues or questions, refer to:
- Main GROMACS skill: `../../SKILL.md`
- Swarm workflow: `../../references/swarm-workflow.md`
- MDP templates: `../../references/mdp-templates.md`
