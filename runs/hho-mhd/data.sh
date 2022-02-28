#!/bin/bash
#

# Test case
u_id=1 #1,2 compatible with Dirichlet, 2 compatible with Hodge, 0 = zero solution, 3 will not work properly
b_id=3 #1,2 compatible with Dirichlet, 2,3 compatible with Hodge, 0 = zero solution
p_id=2 #0 = zero solution, 1 = linear, 2 = sin-sin-sin
u_bc='D'
b_bc='H'
visc=1E0
diff=1E0
tol=1E-6

# Plot File
plotfile="plot"

# Order
k=1
l=1
use_threads="true"

mesh[1]="Voro-small-0/RF_fmt/voro-2"
mesh[2]="Voro-small-0/RF_fmt/voro-4"
mesh[3]="Voro-small-0/RF_fmt/voro-6"
#mesh[4]="Voro-small-0/RF_fmt/voro-8"
#mesh[5]="Voro-small-0/RF_fmt/voro-10"
#mesh[1]="Voro-small-0/RF_fmt/voro-12"

#mesh[1]="Voro-small-1/RF_fmt/voro.1"
#mesh[1]="Voro-small-1/RF_fmt/voro.2"
#mesh[2]="Voro-small-1/RF_fmt/voro.3"
#mesh[3]="Voro-small-1/RF_fmt/voro.4"
#mesh[4]="Voro-small-1/RF_fmt/voro.5"
#mesh[5]="Voro-small-1/RF_fmt/voro.6"
#mesh[6]="Voro-small-1/RF_fmt/voro.7"

#mesh[1]="Voro-small-2/RF_fmt/voro.1"
#mesh[1]="Voro-small-2/RF_fmt/voro.2"
#mesh[2]="Voro-small-2/RF_fmt/voro.3"
#mesh[3]="Voro-small-2/RF_fmt/voro.4"
#mesh[4]="Voro-small-2/RF_fmt/voro.5"
#mesh[5]="Voro-small-2/RF_fmt/voro.6"
#mesh[6]="Voro-small-2/RF_fmt/voro.7"
#mesh[7]="Voro-small-2/RF_fmt/voro.8"

#mesh[1]="Voro-Tets-1/RF_fmt/voro.1"
#mesh[2]="Voro-Tets-1/RF_fmt/voro.2"
#mesh[3]="Voro-Tets-1/RF_fmt/voro.3"
#mesh[4]="Voro-Tets-1/RF_fmt/voro.4"

#mesh[1]="Tetgen-Cube-0/RF_fmt/cube.1"
#mesh[2]="Tetgen-Cube-0/RF_fmt/cube.2"
#mesh[3]="Tetgen-Cube-0/RF_fmt/cube.3"
#mesh[4]="Tetgen-Cube-0/RF_fmt/cube.4"
#mesh[5]="Tetgen-Cube-0/RF_fmt/cube.5"
#mesh[6]="Tetgen-Cube-0/RF_fmt/cube.6"
#mesh[7]="Tetgen-Cube-0/RF_fmt/cube.7"

#mesh[1]="Cubic-Cells/RF_fmt/gcube_1x1x1"
#mesh[1]="Cubic-Cells/RF_fmt/gcube_2x2x2"
#mesh[2]="Cubic-Cells/RF_fmt/gcube_4x4x4"
#mesh[3]="Cubic-Cells/RF_fmt/gcube_8x8x8"
#mesh[3]="Cubic-Cells/RF_fmt/gcube_16x16x16"
#mesh[5]="Cubic-Cells/RF_fmt/gcube_32x32x32"
#mesh[6]="Cubic-Cells/RF_fmt/gcube_64x64x64"





