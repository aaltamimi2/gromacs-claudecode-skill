#!/bin/bash
# GROMACS Wrapper
# Automatically detects and uses gmx_mpi (preferred) or gmx (fallback)

# Detect available GROMACS commands (prioritize gmx_mpi)
if command -v gmx_mpi &>/dev/null; then
    GMX_CMD="gmx_mpi"
elif command -v gmx &>/dev/null; then
    GMX_CMD="gmx"
else
    echo "Error: No GROMACS installation found (gmx_mpi or gmx)" >&2
    exit 1
fi

# Execute the command
exec $GMX_CMD "$@"
