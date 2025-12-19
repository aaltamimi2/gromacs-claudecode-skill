#!/usr/bin/env python3
"""
Convert Gaussian .gjf file to ORCA input format
Handles PDB-style atom names and large coordinates
"""

import sys
import re

def parse_gaussian_gjf(gjf_file):
    """Extract coordinates and settings from Gaussian .gjf file"""
    
    with open(gjf_file, 'r') as f:
        lines = f.readlines()
    
    # Find route section (starts with #)
    route = ""
    for i, line in enumerate(lines):
        if line.strip().startswith('#'):
            route = line.strip()
            break
    
    # Extract charge and multiplicity
    charge_mult_idx = None
    for i, line in enumerate(lines):
        if line.strip() and re.match(r'^[+-]?\d+\s+\d+', line.strip()):
            charge_mult_idx = i
            charge, mult = map(int, line.strip().split()[:2])
            break
    
    if charge_mult_idx is None:
        raise ValueError("Could not find charge and multiplicity")
    
    # Extract coordinates (everything after charge/mult until blank line or connectivity)
    coords = []
    for i in range(charge_mult_idx + 1, len(lines)):
        line = lines[i].strip()
        
        # Stop at blank line or connectivity section
        if not line or line[0].isdigit():
            break
        
        # Parse atom line: Element(PDBName=...) X Y Z L
        # or simple: Element X Y Z
        parts = line.split()
        if len(parts) >= 4:
            # Extract element symbol (before any parenthesis)
            element = parts[0].split('(')[0]
            
            # Get coordinates (last 4 items are: X, Y, Z, L or just X, Y, Z)
            try:
                x = float(parts[-4])
                y = float(parts[-3])
                z = float(parts[-2])
                coords.append((element, x, y, z))
            except:
                # Try without L flag
                x = float(parts[-3])
                y = float(parts[-2])
                z = float(parts[-1])
                coords.append((element, x, y, z))
    
    return {
        'route': route,
        'charge': charge,
        'mult': mult,
        'coords': coords
    }

def write_orca_input(data, output_file, level='HF', basis='3-21G', 
                     nprocs=56, mem_per_core=500):
    """Write ORCA input file"""
    
    with open(output_file, 'w') as f:
        # Header and keywords
        f.write(f"# Converted from Gaussian .gjf\n")
        f.write(f"# Using {level}/{basis}\n\n")
        
        # ORCA keywords (matching HF/3-21G from Gaussian)
        if level.upper() == 'HF':
            f.write(f"! {level} {basis} TightSCF Opt\n\n")
        else:
            f.write(f"! {level} {basis} TightSCF TightOpt Grid5 RI-J\n\n")
        
        # Parallel settings
        f.write(f"%pal nprocs {nprocs} end\n\n")
        f.write(f"%maxcore {mem_per_core}\n\n")
        
        # SCF settings for large molecules
        f.write("%scf\n")
        f.write("  MaxIter 500\n")
        f.write("  DIISMaxEq 10\n")
        f.write("end\n\n")
        
        # Coordinates
        f.write(f"* xyz {data['charge']} {data['mult']}\n")
        for element, x, y, z in data['coords']:
            f.write(f"  {element:2s}  {x:14.8f}  {y:14.8f}  {z:14.8f}\n")
        f.write("*\n")

def write_orca_xyz(data, xyz_file):
    """Write XYZ file for later use"""
    
    with open(xyz_file, 'w') as f:
        f.write(f"{len(data['coords'])}\n")
        f.write(f"Converted from Gaussian .gjf, charge={data['charge']} mult={data['mult']}\n")
        for element, x, y, z in data['coords']:
            f.write(f"{element:2s}  {x:14.8f}  {y:14.8f}  {z:14.8f}\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python gjf_to_orca.py input.gjf [output_base_name]")
        print("  Creates: output_base_name.inp and output_base_name.xyz")
        sys.exit(1)
    
    gjf_file = sys.argv[1]
    
    # Output base name
    if len(sys.argv) > 2:
        base_name = sys.argv[2]
    else:
        base_name = gjf_file.replace('.gjf', '')
    
    print(f"Reading Gaussian file: {gjf_file}")
    
    try:
        # Parse Gaussian file
        data = parse_gaussian_gjf(gjf_file)
        
        print(f"  Found {len(data['coords'])} atoms")
        print(f"  Charge: {data['charge']}, Multiplicity: {data['mult']}")
        
        # Write ORCA input
        orca_file = f"{base_name}.inp"
        xyz_file = f"{base_name}.xyz"
        
        write_orca_input(data, orca_file)
        write_orca_xyz(data, xyz_file)
        
        print(f"\nCreated ORCA files:")
        print(f"  Input: {orca_file}")
        print(f"  XYZ:   {xyz_file}")
        print(f"\nTo run: orca {orca_file} > {base_name}.out")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
