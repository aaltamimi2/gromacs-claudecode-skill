# GROMACS Analysis Workflow Templates

## Template 1: Standard Protein Analysis

```bash
#!/bin/bash
# Standard protein MD analysis pipeline

TRAJ="md.xtc"
TPR="md.tpr"
PREFIX="analysis"

# 1. RMSD (backbone)
echo "4\n4\n" | gmx rms -f $TRAJ -s $TPR -o ${PREFIX}_rmsd_backbone.xvg -tu ns

# 2. RMSD (C-alpha)
echo "3\n3\n" | gmx rms -f $TRAJ -s $TPR -o ${PREFIX}_rmsd_ca.xvg -tu ns

# 3. RMSF (C-alpha)
echo "3\n" | gmx rmsf -f $TRAJ -s $TPR -o ${PREFIX}_rmsf_ca.xvg -res

# 4. Radius of gyration
echo "1\n" | gmx gyrate -f $TRAJ -s $TPR -o ${PREFIX}_gyrate.xvg

# 5. SASA
echo "1\n" | gmx sasa -f $TRAJ -s $TPR -o ${PREFIX}_sasa.xvg -tu ns

# 6. Secondary structure
echo "1\n" | gmx do_dssp -f $TRAJ -s $TPR -o ${PREFIX}_ss.xpm -sc ${PREFIX}_scount.xvg

# 7. H-bonds
echo "1\n1\n" | gmx hbond -f $TRAJ -s $TPR -num ${PREFIX}_hbnum.xvg

echo "Analysis complete! Files: ${PREFIX}_*.xvg"
```

## Template 2: Membrane Protein Analysis

```bash
#!/bin/bash
# Membrane protein specific analysis

TRAJ="md.xtc"
TPR="md.tpr"

# 1. Membrane thickness
echo "Lipids\n" | gmx density -f $TRAJ -s $TPR -o membrane_density.xvg -d Z -sl 100

# 2. Protein tilt relative to membrane normal
gmx gangle -f $TRAJ -s $TPR -g1 vector -group1 'Protein' -g2 vector -group2 'Lipids' -oav tilt.xvg

# 3. Lipid order parameters
gmx order -f $TRAJ -s $TPR -n index.ndx -od order.xvg

# 4. Protein-lipid contacts
echo "Protein\nLipids\n" | gmx mindist -f $TRAJ -s $TPR -od protein_lipid_dist.xvg
```

## Template 3: Ligand Binding Analysis

```bash
#!/bin/bash
# Protein-ligand specific analysis

TRAJ="md.xtc"
TPR="md.tpr"
LIG="LIG"  # Ligand residue name

# 1. Ligand RMSD
echo "resname $LIG\nresname $LIG\n" | gmx rms -f $TRAJ -s $TPR -o ligand_rmsd.xvg

# 2. Distance from binding site
# (Replace 'resid 45' with your binding site residue)
gmx distance -f $TRAJ -s $TPR -select "com of resname $LIG" "com of resid 45" -oall ligand_distance.xvg

# 3. H-bonds between protein and ligand
echo "Protein\nresname $LIG\n" | gmx hbond -f $TRAJ -s $TPR -num protein_lig_hbonds.xvg

# 4. Ligand SASA (burial)
echo "resname $LIG\n" | gmx sasa -f $TRAJ -s $TPR -o ligand_sasa.xvg
```

## Template 4: Quick Quality Check

```bash
#!/bin/bash
# Quick simulation quality check

TPR="md.tpr"
EDR="md.edr"
TRAJ="md.xtc"

echo "=== Quick MD Quality Check ==="

# Energy conservation
echo "Total-Energy\n0\n" | gmx energy -f $EDR -o energy.xvg
echo "✓ Check energy.xvg for drift"

# Temperature
echo "Temperature\n0\n" | gmx energy -f $EDR -o temperature.xvg
TEMP=$(tail -1 temperature.xvg | awk '{print $2}')
echo "✓ Average temperature: $TEMP K"

# Pressure (if NPT)
echo "Pressure\n0\n" | gmx energy -f $EDR -o pressure.xvg 2>/dev/null && echo "✓ Pressure data extracted"

# Density (if NPT)
echo "Density\n0\n" | gmx energy -f $EDR -o density.xvg 2>/dev/null
DENS=$(tail -1 density.xvg | awk '{print $2}')
echo "✓ Average density: $DENS kg/m³"

# Trajectory integrity
gmx check -f $TRAJ
echo "✓ Trajectory integrity checked"

echo "=== Quality check complete ==="
```

## Using with Claude Code

Save these as templates and ask Claude:
- "Run the standard protein analysis workflow on my trajectory"
- "Adapt the membrane analysis template for my GPCR simulation"
- "Create a custom analysis combining RMSD, RMSF, and H-bonds"
