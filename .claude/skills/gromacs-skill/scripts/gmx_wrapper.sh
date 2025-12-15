#!/bin/bash
# GROMACS MPI Wrapper
# Automatically detects whether to use mpirun for gmx_mpi commands

# Detect available GROMACS commands
if command -v gmx &>/dev/null; then
    GMX_CMD="gmx"
elif command -v gmx_mpi &>/dev/null; then
    # For single-process runs, use mpirun -np 1
    GMX_CMD="mpirun -np 1 gmx_mpi"
else
    echo "Error: No GROMACS installation found (gmx or gmx_mpi)" >&2
    exit 1
fi

# Execute the command
$GMX_CMD "$@"
