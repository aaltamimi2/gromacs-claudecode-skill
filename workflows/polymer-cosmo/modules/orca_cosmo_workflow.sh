#!/bin/bash
###############################################################################
# Complete ORCA COSMO-RS Workflow (matching Gaussian BVP86/TZVP workflow)
# 
# This replicates the Gaussian workflow:
#   1. Gas phase optimization (BP86/TZVP)
#   2. CPCM optimization in water (BP86/TZVP)
#   3. COSMO-RS single point (infinite dielectric)
#
# Usage: ./orca_cosmo_workflow.sh molecule_name.xyz [nprocs]
###############################################################################

set -e  # Exit on error

# Input validation
if [ $# -eq 0 ]; then
    echo "Usage: $0 molecule.xyz [nprocs]"
    echo ""
    echo "Example: $0 conf_2300.xyz 56"
    echo ""
    echo "This will create:"
    echo "  - molecule_gas.inp/out     (gas phase optimization)"
    echo "  - molecule_cpcm.inp/out    (CPCM optimization in water)"
    echo "  - molecule_cosmo.inp/out   (COSMO-RS single point)"
    echo "  - molecule_cosmo.cosmo     (COSMO file for COSMOtherm)"
    exit 1
fi

XYZ_FILE=$1
NPROCS=${2:-56}  # Default to 56 cores

# Extract base name (remove .xyz extension)
BASE_NAME=$(basename "$XYZ_FILE" .xyz)

# Check if XYZ file exists
if [ ! -f "$XYZ_FILE" ]; then
    echo "Error: File $XYZ_FILE not found!"
    exit 1
fi

# Memory per core (MB) - adjust based on your system
MEM_PER_CORE=500

echo "###############################################################################"
echo "ORCA COSMO-RS Workflow"
echo "###############################################################################"
echo "Molecule:        $BASE_NAME"
echo "Input XYZ:       $XYZ_FILE"
echo "Cores:           $NPROCS"
echo "Memory/core:     ${MEM_PER_CORE} MB"
echo "###############################################################################"
echo ""

# Extract charge and multiplicity from XYZ comment line
# Format: "comment line, charge=X mult=Y" or just assume 0 1
CHARGE=0
MULT=1
if grep -q "charge=" "$XYZ_FILE"; then
    CHARGE=$(grep -oP 'charge=\K\d+' "$XYZ_FILE" | head -1)
    MULT=$(grep -oP 'mult=\K\d+' "$XYZ_FILE" | head -1)
fi

###############################################################################
# Step 1: Gas Phase Optimization
###############################################################################
echo "============================"
echo "Step 1: Gas Phase Optimization"
echo "============================"
echo "Method: BP86/def2-TZVP"
echo "Started: $(date)"
echo ""

GAS_INP="${BASE_NAME}_gas.inp"
GAS_OUT="${BASE_NAME}_gas.out"
GAS_XYZ="${BASE_NAME}_gas.xyz"

cat > "$GAS_INP" << EOF
# Gas phase optimization with BP86/def2-TZVP
# Equivalent to Gaussian: # bvp86/tzvp/dga1 scf=tight opt=tight

! BP86 def2-TZVP def2/J TightSCF TightOpt Grid5
! RI-J

%pal nprocs $NPROCS end
%maxcore $MEM_PER_CORE

%scf
  MaxIter 500
  DIISMaxEq 10
end

%geom
  MaxIter 500
  TolE 5e-6      # Tight energy convergence
  TolRMSG 1e-4   # Tight gradient convergence  
  TolMaxG 3e-4
end

* xyzfile $CHARGE $MULT $XYZ_FILE
EOF

echo "Running gas phase optimization..."
orca "$GAS_INP" > "$GAS_OUT"

if grep -q "ORCA TERMINATED NORMALLY" "$GAS_OUT"; then
    E_GAS=$(grep "FINAL SINGLE POINT ENERGY" "$GAS_OUT" | tail -1 | awk '{print $5}')
    echo "✓ Gas phase optimization complete"
    echo "  Energy: $E_GAS Eh"
    echo "  Optimized geometry: $GAS_XYZ"
else
    echo "✗ Gas phase optimization FAILED"
    echo "  Check $GAS_OUT for errors"
    exit 1
fi

echo ""

###############################################################################
# Step 2: CPCM Optimization in Water
###############################################################################
echo "============================"
echo "Step 2: CPCM Optimization"
echo "============================"
echo "Method: BP86/def2-TZVP with CPCM(Water)"
echo "Started: $(date)"
echo ""

CPCM_INP="${BASE_NAME}_cpcm.inp"
CPCM_OUT="${BASE_NAME}_cpcm.out"
CPCM_XYZ="${BASE_NAME}_cpcm.xyz"

cat > "$CPCM_INP" << EOF
# CPCM optimization in water with BP86/def2-TZVP
# Equivalent to Gaussian: # bvp86/tzvp/dga1 opt=tight scf=tight SCRF=(CPCM)

! BP86 def2-TZVP def2/J TightSCF TightOpt Grid5
! CPCM(Water) RI-J

%pal nprocs $NPROCS end
%maxcore $MEM_PER_CORE

%scf
  MaxIter 500
  DIISMaxEq 10
end

%geom
  MaxIter 500
  TolE 5e-6
  TolRMSG 1e-4
  TolMaxG 3e-4
end

* xyzfile $CHARGE $MULT $GAS_XYZ
EOF

echo "Running CPCM optimization..."
orca "$CPCM_INP" > "$CPCM_OUT"

if grep -q "ORCA TERMINATED NORMALLY" "$CPCM_OUT"; then
    E_CPCM=$(grep "FINAL SINGLE POINT ENERGY" "$CPCM_OUT" | tail -1 | awk '{print $5}')
    echo "✓ CPCM optimization complete"
    echo "  Energy: $E_CPCM Eh"
    echo "  Optimized geometry: $CPCM_XYZ"
else
    echo "✗ CPCM optimization FAILED"
    echo "  Check $CPCM_OUT for errors"
    exit 1
fi

echo ""

###############################################################################
# Step 3: COSMO-RS Single Point (Infinite Dielectric)
###############################################################################
echo "============================"
echo "Step 3: COSMO-RS Single Point"
echo "============================"
echo "Method: BP86/def2-TZVP with COSMO (epsilon=infinity)"
echo "Started: $(date)"
echo ""

COSMO_INP="${BASE_NAME}_cosmo.inp"
COSMO_OUT="${BASE_NAME}_cosmo.out"
COSMO_FILE="${BASE_NAME}_cosmo.cosmo"

cat > "$COSMO_INP" << EOF
# COSMO-RS single point calculation
# Equivalent to Gaussian: # bvp86/tzvp/dga1 geom=checkpoint guess=read scf=tight SCRF=COSMORS

! BP86 def2-TZVP def2/J TightSCF Grid5
! RI-J

%pal nprocs $NPROCS end
%maxcore $MEM_PER_CORE

%scf
  MaxIter 500
  DIISMaxEq 10
end

%cpcm
  epsilon infinity  # Infinite dielectric constant for COSMO-RS
end

* xyzfile $CHARGE $MULT $CPCM_XYZ
EOF

echo "Running COSMO-RS single point..."
orca "$COSMO_INP" > "$COSMO_OUT"

if grep -q "ORCA TERMINATED NORMALLY" "$COSMO_OUT"; then
    E_COSMO=$(grep "FINAL SINGLE POINT ENERGY" "$COSMO_OUT" | tail -1 | awk '{print $5}')
    echo "✓ COSMO-RS single point complete"
    echo "  Energy: $E_COSMO Eh"
    
    # Check for COSMO file
    if [ -f "$COSMO_FILE" ]; then
        echo "  ✓ COSMO file generated: $COSMO_FILE"
    else
        echo "  ✗ WARNING: COSMO file not found!"
        echo "    Expected: $COSMO_FILE"
    fi
else
    echo "✗ COSMO-RS single point FAILED"
    echo "  Check $COSMO_OUT for errors"
    exit 1
fi

echo ""

###############################################################################
# Summary
###############################################################################
echo "###############################################################################"
echo "WORKFLOW COMPLETE"
echo "###############################################################################"
echo ""
echo "Results Summary:"
echo "  Molecule:          $BASE_NAME"
echo "  Gas phase energy:  $E_GAS Eh"
echo "  CPCM energy:       $E_CPCM Eh"
echo "  COSMO energy:      $E_COSMO Eh"
echo ""

# Calculate solvation free energy (if bc is available)
if command -v bc &> /dev/null; then
    # ΔG_solv = (E_CPCM - E_gas) × 627.509 kcal/mol
    DG_SOLV=$(echo "scale=2; ($E_CPCM - $E_GAS) * 627.509" | bc)
    echo "  Solvation ΔG:      $DG_SOLV kcal/mol"
    echo ""
fi

echo "Generated files:"
echo "  Gas phase:         $GAS_INP, $GAS_OUT, $GAS_XYZ"
echo "  CPCM:              $CPCM_INP, $CPCM_OUT, $CPCM_XYZ"
echo "  COSMO-RS:          $COSMO_INP, $COSMO_OUT"
echo "  COSMO file:        $COSMO_FILE"
echo ""
echo "Next step: Use $COSMO_FILE with COSMOtherm for predictions"
echo "###############################################################################"
