#!/bin/bash
#

# Face and cell degrees
k=2
l=1

# Use threads
use_threads="true"

# Boundary conditions (D, N, Mx)
bc="M0"

# Test case
tcsol=1
tcdiff=2

# Solver
solver_type="bicgstab"
#solver_type="ma41"

# Meshes (with type)
#mesh[1]="RF:Tetgen-Cube-0/RF_fmt/cube.1"
#mesh[2]="RF:Tetgen-Cube-0/RF_fmt/cube.2"
#mesh[3]="RF:Tetgen-Cube-0/RF_fmt/cube.3"
#mesh[4]="RF:Tetgen-Cube-0/RF_fmt/cube.4"
#mesh[5]="RF:Tetgen-Cube-0/RF_fmt/cube.5"
#mesh[6]="RF:Tetgen-Cube-0/RF_fmt/cube.6"

#mesh[1]="RF:Tetgen-Cube-1/RF_fmt/cube.2"
#mesh[2]="RF:Tetgen-Cube-1/RF_fmt/cube.3"
#mesh[3]="RF:Tetgen-Cube-1/RF_fmt/cube.4"
#mesh[4]="RF:Tetgen-Cube-1/RF_fmt/cube.5"

#mesh[1]="RF:Random-Hexahedra/RF_fmt/gcube.1"
#mesh[2]="RF:Random-Hexahedra/RF_fmt/gcube.2"
#mesh[3]="RF:Random-Hexahedra/RF_fmt/gcube.3"

#mesh[1]="RF:Prysmatic-Cells/RF_fmt/vprism_4x4x4"
#mesh[2]="RF:Prysmatic-Cells/RF_fmt/vprism_8x8x8"
#mesh[3]="RF:Prysmatic-Cells/RF_fmt/vprism_16x16x16"

#mesh[1]="RF:Cubic-Cells/RF_fmt/gcube_4x4x4"
#mesh[2]="RF:Cubic-Cells/RF_fmt/gcube_8x8x8"
#mesh[3]="RF:Cubic-Cells/RF_fmt/gcube_16x16x16"
#mesh[4]="RF:Cubic-Cells/RF_fmt/gcube_32x32x32"

#mesh[1]="RF:Voro-Tets-2/RF_fmt/voro.1"
#mesh[2]="RF:Voro-Tets-2/RF_fmt/voro.2"
#mesh[3]="RF:Voro-Tets-2/RF_fmt/voro.3"
#mesh[4]="RF:Voro-Tets-2/RF_fmt/voro.4"
#mesh[5]="RF:Voro-Tets-2/RF_fmt/voro.5"

#mesh[1]="RF:Voro-Tets-1/RF_fmt/voro.1"
#mesh[2]="RF:Voro-Tets-1/RF_fmt/voro.2"
#mesh[3]="RF:Voro-Tets-1/RF_fmt/voro.3"
#mesh[4]="RF:Voro-Tets-1/RF_fmt/voro.4"
#mesh[5]="RF:Voro-Tets-1/RF_fmt/voro.5"
#mesh[6]="RF:Voro-Tets-1/RF_fmt/voro.6"

mesh[1]="RF:Voro-small-0/RF_fmt/voro-2"
mesh[2]="RF:Voro-small-0/RF_fmt/voro-4"
mesh[3]="RF:Voro-small-0/RF_fmt/voro-6"
#mesh[4]="RF:Voro-small-0/RF_fmt/voro-8"
#mesh[5]="RF:Voro-small-0/RF_fmt/voro-10"

#mesh[1]="RF:Voro-small-1/RF_fmt/voro.3"
#mesh[2]="RF:Voro-small-1/RF_fmt/voro.4"
#mesh[3]="RF:Voro-small-1/RF_fmt/voro.5"
#mesh[4]="RF:Voro-small-1/RF_fmt/voro.6"

#mesh[1]="RF:Voro-small-2/RF_fmt/voro.5"
#mesh[2]="RF:Voro-small-2/RF_fmt/voro.6"
#mesh[3]="RF:Voro-small-2/RF_fmt/voro.7"
#mesh[4]="RF:Voro-small-2/RF_fmt/voro.8"

#mesh[1]="RF:Voro-Cube-0/RF_fmt/voro.5"
#mesh[2]="RF:Voro-Cube-0/RF_fmt/voro.10"
#mesh[3]="RF:Voro-Cube-0/RF_fmt/voro.20"
#mesh[4]="RF:Voro-Cube-0/RF_fmt/voro.30"

