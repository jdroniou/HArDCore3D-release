#!/bin/bash
#

# Degree
k=1

# Stabilization parameters
stab_par_stokes=3
stab_par_darcy=.3

# Test case and scalings of viscosity, permeability
tcsol=2
scaling_viscosity=1
scaling_permeabilityinv=1

# Mesh family
 mesh_family=Tetgen-Cube-0
# mesh_family=Cubic-Cells
# mesh_family=Tetgen-Cube-1
# mesh_family=Random-Hexahedra
# mesh_family=Prysmatic-Cells
# mesh_family=Voro-Tets-1
# mesh_family=Voro-Tets-2
#mesh_family=Voro-small-0
# mesh_family=Voro-small-1
# mesh_family=Voro-small-2
# mesh_family=Voro-Cube-0
#mesh_family=CubeInCube-Tet
#mesh_family=CubeCavityWedge-hexa
#mesh_family=CubeCavityWedge-tet-hexa

case ${mesh_family} in
    Tetgen-Cube-0)
	mesh[1]="RF:Tetgen-Cube-0/RF_fmt/cube.1"
	mesh[2]="RF:Tetgen-Cube-0/RF_fmt/cube.2"
	mesh[3]="RF:Tetgen-Cube-0/RF_fmt/cube.3"
	mesh[4]="RF:Tetgen-Cube-0/RF_fmt/cube.4"
	mesh[5]="RF:Tetgen-Cube-0/RF_fmt/cube.5"
	mesh[6]="RF:Tetgen-Cube-0/RF_fmt/cube.6"
	;;
    Tetgen-Cube-1)
	mesh[1]="RF:Tetgen-Cube-1/RF_fmt/cube.2"
	mesh[2]="RF:Tetgen-Cube-1/RF_fmt/cube.3"
	mesh[3]="RF:Tetgen-Cube-1/RF_fmt/cube.4"
	mesh[4]="RF:Tetgen-Cube-1/RF_fmt/cube.5"
	;;
    Random-Hexahedra)
	mesh[1]="RF:Random-Hexahedra/RF_fmt/gcube.1"
	mesh[2]="RF:Random-Hexahedra/RF_fmt/gcube.2"
#	mesh[3]="RF:Random-Hexahedra/RF_fmt/gcube.3"
	;;
    Prysmatic-Cells)
	mesh[1]="RF:Prysmatic-Cells/RF_fmt/vprism_4x4x4"
	mesh[2]="RF:Prysmatic-Cells/RF_fmt/vprism_8x8x8"
	mesh[3]="RF:Prysmatic-Cells/RF_fmt/vprism_16x16x16"
	;;
    Cubic-Cells)
	mesh[1]="RF:Cubic-Cells/RF_fmt/gcube_2x2x2"
	mesh[2]="RF:Cubic-Cells/RF_fmt/gcube_4x4x4"
	mesh[3]="RF:Cubic-Cells/RF_fmt/gcube_8x8x8"
#	mesh[4]="RF:Cubic-Cells/RF_fmt/gcube_16x16x16"
	# mesh[5]="RF:Cubic-Cells/RF_fmt/gcube_32x32x32"
	;;
    Voro-Tets-2)
	mesh[1]="RF:Voro-Tets-2/RF_fmt/voro.1"
	mesh[2]="RF:Voro-Tets-2/RF_fmt/voro.2"
	mesh[3]="RF:Voro-Tets-2/RF_fmt/voro.3"
	mesh[4]="RF:Voro-Tets-2/RF_fmt/voro.4"
	mesh[5]="RF:Voro-Tets-2/RF_fmt/voro.5"
	;;
    Voro-Tets-1)
	mesh[1]="RF:Voro-Tets-1/RF_fmt/voro.1"
	mesh[2]="RF:Voro-Tets-1/RF_fmt/voro.2"
	mesh[3]="RF:Voro-Tets-1/RF_fmt/voro.3"
	mesh[4]="RF:Voro-Tets-1/RF_fmt/voro.4"
	mesh[5]="RF:Voro-Tets-1/RF_fmt/voro.5"
	mesh[6]="RF:Voro-Tets-1/RF_fmt/voro.6"
	;;
    Voro-small-0)
	mesh[1]="RF:Voro-small-0/RF_fmt/voro-2"
	mesh[2]="RF:Voro-small-0/RF_fmt/voro-4"
  mesh[3]="RF:Voro-small-0/RF_fmt/voro-6"
#	mesh[4]="RF:Voro-small-0/RF_fmt/voro-8"
#	mesh[5]="RF:Voro-small-0/RF_fmt/voro-10"
#	mesh[6]="RF:Voro-small-0/RF_fmt/voro-12"
	# mesh[7]="RF:Voro-small-0/RF_fmt/voro-14"
	# mesh[8]="RF:Voro-small-0/RF_fmt/voro-16"
	;;
    Voro-small-1)
	mesh[1]="RF:Voro-small-1/RF_fmt/voro.3"
	mesh[2]="RF:Voro-small-1/RF_fmt/voro.4"
	mesh[3]="RF:Voro-small-1/RF_fmt/voro.5"
	mesh[4]="RF:Voro-small-1/RF_fmt/voro.6"
	;;
    Voro-small-2)
	mesh[1]="RF:Voro-small-2/RF_fmt/voro.5"
	mesh[2]="RF:Voro-small-2/RF_fmt/voro.6"
	mesh[3]="RF:Voro-small-2/RF_fmt/voro.7"
	mesh[4]="RF:Voro-small-2/RF_fmt/voro.8"
	;;
    Voro-Cube-0)
	mesh[1]="RF:Voro-Cube-0/RF_fmt/voro.5"
	mesh[2]="RF:Voro-Cube-0/RF_fmt/voro.10"
	mesh[3]="RF:Voro-Cube-0/RF_fmt/voro.20"
	mesh[4]="RF:Voro-Cube-0/RF_fmt/voro.30"
	;;
	  CubeInCube-Tet)
  mesh[1]="RF:CubeInCube-Tet/RF_fmt/CIC-tet-0"
  mesh[2]="RF:CubeInCube-Tet/RF_fmt/CIC-tet-1"
  mesh[3]="RF:CubeInCube-Tet/RF_fmt/CIC-tet-2"
  ;;
    CubeCavityWedge-tet-hexa)
  mesh[1]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/tet-hexa.1"
  mesh[2]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/tet-hexa.2"
  mesh[3]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/tet-hexa.3"
#  mesh[4]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/tet-hexa.4"
#  mesh[5]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/tet-hexa.5"
  ;;
esac

