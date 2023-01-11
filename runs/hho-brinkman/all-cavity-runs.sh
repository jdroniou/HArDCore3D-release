#!/bin/bash
#
# Execute hho-brinkman on series of mesh families and degrees, for the cavity test case (solution=6)
#

# Options
if [[ $1 == "help" ]]; then
    echo -e "\nExecute tests using parameters in data.sh"
    exit;
fi;

# Load data: only for scaling and stabilisation, really
. data.sh
tcsol=6

# File for times
timesfile="times.dat"

# Directories
executable_name="hho-brinkman"

origin=$(pwd)
if [ ! -f ../directories.sh ]; then
  echo "directories.sh does not exist. Please read the README.txt in the parent folder."
  exit
fi
. ../directories.sh

# Mesh families and maximum degree
#mesh_families="Voro-small-0 Tetgen-Cube-0"
#mesh_families="Tetgen-Cube-0"
#mesh_families="CubeCavityWedge-hexa"
mesh_families="CubeCavityWedge-tet-hexa"

for mesh_family in ${mesh_families}
do

    unset mesh;
    
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
	          mesh[3]="RF:Random-Hexahedra/RF_fmt/gcube.3"
	          mesh[4]="RF:Random-Hexahedra/RF_fmt/gcube.4"
#	          mesh[5]="RF:Random-Hexahedra/RF_fmt/gcube.5"
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
            mesh[4]="RF:Cubic-Cells/RF_fmt/gcube_16x16x16"
	          # mesh[5]="RF:Cubic-Cells/RF_fmt/gcube_32x32x32"
	          ;;
        Cubic-Cells-coarse)
            mesh[1]="RF:Cubic-Cells/RF_fmt/gcube_16x16x16.coarse.5"
            mesh[2]="RF:Cubic-Cells/RF_fmt/gcube_16x16x16.coarse.3"
            mesh[3]="RF:Cubic-Cells/RF_fmt/gcube_16x16x16.coarse.2"
            mesh[4]="RF:Cubic-Cells/RF_fmt/gcube_16x16x16.coarse.1"
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
	          mesh[4]="RF:Voro-small-0/RF_fmt/voro-8"
	          mesh[5]="RF:Voro-small-0/RF_fmt/voro-10"
	          mesh[6]="RF:Voro-small-0/RF_fmt/voro-12"
	          mesh[7]="RF:Voro-small-0/RF_fmt/voro-14"
	          mesh[8]="RF:Voro-small-0/RF_fmt/voro-16"
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
        CubeCavityWedge-hexa)
	          mesh[1]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/hexa.1"
	          mesh[2]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/hexa.2"
	          mesh[3]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/hexa.3"
	          mesh[4]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/hexa.4"
	          mesh[5]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/hexa.5"
            ;;
        CubeCavityWedge-tet-hexa)
	          mesh[1]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/tet-hexa.1"
	          mesh[2]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/tet-hexa.2"
	          mesh[3]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/tet-hexa.3"
	          mesh[4]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/tet-hexa.4"
	          mesh[5]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/tet-hexa.5"
            ;;
    esac

    for k in {0..2}
    do
      
      echo "Degree                 : $k"
      echo "Mesh family            : $mesh_family"
      echo "Stabilization parameter: $stab_par"
      echo "Test case: solution    : $tcsol"
      echo "Output directory       : $outdir"

      outsubdir=${outdir}/${mesh_family}_k${k}
      if [ ! -d $outsubdir ]; then
          mkdir -p $outsubdir
      else
          \rm -r $outsubdir/*
      fi

      ###
      # EXECUTE FOR EACH MESH
      nbmesh=${#mesh[@]}
      for i in `seq 1 $nbmesh`; 
      do
          meshfile=$meshdir"/"$(echo ${mesh[$i]} | cut -d ':' -f 2)
          echo -e "------------------------------------------------------------------------------"
          echo -e "Mesh $i out of $nbmesh: $meshfile"
          echo -e "Output directory: $outsubdir"
          echo -e "------------------------------------------------------------------------------"
          # Execute code
          if($executable -m $meshfile -k $k -s $tcsol -x $stab_par_stokes -y $stab_par_darcy -v $scaling_viscosity -p $scaling_permeabilityinv); then
            # Move outputs
            mv results.txt $outsubdir/results-$i.txt
          fi
      done

      # CREATE FILE FOR FLUXES (not really useful actually, will need to run separate script once all the tests are in)
#      fluxfile=flux_convergence.dat;
#      refFlux=-3e-7; # To be taken from smallest h/largest k
#      echo -e "Deg MeshSize NbCells NbFaces SizeSystem Flux Rate" > $outsubdir/$fluxfile
#          for i in `seq 1 $nbmesh`; 
#          do
#              Degree=$(awk '/Degree:/ {print $NF}' $outsubdir/results-$i.txt)
#              MeshSize=$(awk '/MeshSize:/ {print $NF}' $outsubdir/results-$i.txt)
#              NbCells=$(awk '/NbCells:/ {print $NF}' $outsubdir/results-$i.txt)
#              NbFaces=$(awk '/NbFaces:/ {print $NF}' $outsubdir/results-$i.txt)
#              SizeSystem=$(awk '/SizeSystem:/ {print $NF}' $outsubdir/results-$i.txt)
#              Flux=$(awk '/FluxCavity:/ {print $NF}' $outsubdir/results-$i.txt)
#              FluxDiff=$(perl -E "say $Flux - $refFlux");
#              if(($i > 1)); then
#                imo=$(perl -E "say $i - 1")
#                OldMeshSize=$(awk '/MeshSize:/ {print $NF}' $outsubdir/results-$imo.txt)
#                OldFlux=$(awk '/FluxCavity:/ {print $NF}' $outsubdir/results-$imo.txt)
#                OldFluxDiff=$(perl -E "say $OldFlux - $refFlux");
#                FluxRate=$(perl -E "say sprintf(\"%.2f\", log($OldFluxDiff/$FluxDiff)/log($OldMeshSize/$MeshSize))")
#                echo -e "$Degree $MeshSize $NbCells $NbFaces $SizeSystem $FluxDiff $FluxRate" >> $outsubdir/$fluxfile
#              else
#            	  echo -e "$Degree $MeshSize $NbCells $NbFaces $SizeSystem $FluxDiff --" >> $outsubdir/$fluxfile
#              fi
#          done;
#          column -t < $outsubdir/$fluxfile | tee $outsubdir/$fluxfile


    done # for degree

done # for mesh_family

