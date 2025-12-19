#!/usr/bin/env python3
"""
Convert PDB file to XYZ format for ORCA
Handles both standard PDB and Gaussian-style PDB coordinates
"""

import sys
import argparse

def parse_pdb(pdb_file):
    """Extract coordinates from PDB file"""
    
    atoms = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            # Standard PDB ATOM/HETATM lines
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # PDB format:
                # ATOM      1  N   MET A   1      20.154  29.699   5.276  1.00 49.05           N
                # Columns: 0-6(ATOM), 6-11(serial), 12-16(name), 17(altLoc), 17-20(resName),
                #          21(chainID), 22-26(resSeq), 27(iCode), 30-38(x), 38-46(y), 46-54(z),
                #          54-60(occupancy), 60-66(tempFactor), 76-78(element)
                
                try:
                    element = line[76:78].strip()
                    if not element:
                        # Try atom name if element not specified
                        atom_name = line[12:16].strip()
                        element = ''.join([c for c in atom_name if c.isalpha()])[:2]
                    
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    
                    atoms.append((element, x, y, z))
                    
                except (ValueError, IndexError):
                    # Try alternate parsing for non-standard PDB
                    parts = line.split()
                    if len(parts) >= 6:
                        element = parts[2] if parts[2].isalpha() else parts[1]
                        element = ''.join([c for c in element if c.isalpha()])
                        x, y, z = float(parts[-6]), float(parts[-5]), float(parts[-4])
                        atoms.append((element, x, y, z))
    
    return atoms

def write_xyz(atoms, xyz_file, title="Converted from PDB", charge=0, mult=1):
    """Write XYZ file with charge/mult in comment"""
    
    with open(xyz_file, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{title}, charge={charge} mult={mult}\n")
        for element, x, y, z in atoms:
            f.write(f"{element:2s}  {x:14.8f}  {y:14.8f}  {z:14.8f}\n")

def main():
    parser = argparse.ArgumentParser(
        description='Convert PDB file to XYZ format for ORCA',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic conversion
  python pdb_to_xyz.py molecule.pdb
  
  # Specify charge and multiplicity
  python pdb_to_xyz.py molecule.pdb -c 1 -m 2
  
  # Custom output name
  python pdb_to_xyz.py molecule.pdb -o output.xyz
        """
    )
    
    parser.add_argument('pdb_file', help='Input PDB file')
    parser.add_argument('-o', '--output', help='Output XYZ file (default: input_name.xyz)')
    parser.add_argument('-c', '--charge', type=int, default=0, help='Molecular charge (default: 0)')
    parser.add_argument('-m', '--mult', type=int, default=1, help='Spin multiplicity (default: 1)')
    parser.add_argument('-t', '--title', default=None, help='Title for XYZ file')
    
    args = parser.parse_args()
    
    # Determine output filename
    if args.output:
        xyz_file = args.output
    else:
        xyz_file = args.pdb_file.replace('.pdb', '.xyz')
        if xyz_file == args.pdb_file:
            xyz_file = args.pdb_file + '.xyz'
    
    # Title
    if args.title:
        title = args.title
    else:
        title = f"Converted from {args.pdb_file}"
    
    print(f"Reading PDB file: {args.pdb_file}")
    
    try:
        # Parse PDB
        atoms = parse_pdb(args.pdb_file)
        
        if not atoms:
            print("Error: No atoms found in PDB file!")
            sys.exit(1)
        
        print(f"  Found {len(atoms)} atoms")
        
        # Element statistics
        from collections import Counter
        elements = Counter([atom[0] for atom in atoms])
        print("  Composition:", ', '.join([f"{elem}: {count}" for elem, count in sorted(elements.items())]))
        
        # Write XYZ
        write_xyz(atoms, xyz_file, title, args.charge, args.mult)
        
        print(f"\nCreated XYZ file: {xyz_file}")
        print(f"  Charge: {args.charge}")
        print(f"  Multiplicity: {args.mult}")
        print(f"\nNext step: ./orca_cosmo_workflow.sh {xyz_file}")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
