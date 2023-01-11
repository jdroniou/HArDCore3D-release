#!/bin/bash
#
# Execute hho-brinkman on series of mesh families and degrees, and calculate outputs
#

# Options
if [[ $1 == "help" ]]; then
    echo -e "\nExecute tests using parameters in data.sh"
    exit;
fi;

# Load data: meshes and degrees will be replaced, only the other parameters are kept
. data.sh

# File for times
timesfile="times.dat"

# Directories
executable_name="hho-brinkman"
#make $executable_name

origin=$(pwd)
if [ ! -f ../directories.sh ]; then
  echo "directories.sh does not exist. Please read the README.txt in the parent folder."
  exit
fi
. ../directories.sh

# Mesh families and maximum degree
#mesh_families="Voro-small-0 Tetgen-Cube-0"
mesh_families="Random-Hexahedra"
#mesh_families="CubeCavityWedge-tet-hexa"

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
        CubeCavityWedge-tet-hexa)
	          mesh[1]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/tet-hexa.1"
	          mesh[2]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/tet-hexa.2"
	          mesh[3]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/tet-hexa.3"
	          mesh[4]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/tet-hexa.4"
	          mesh[5]="RF:ComplexGeometries/CubeCavityWedge/RF_fmt/tet-hexa.5"
            ;;
    esac

    for k in {1..2}
    do
      
      for scaling_viscosity in {0..1}
      do
      
       for scaling_permeabilityinv in  {0..1}
        do
   
         test=$(($scaling_viscosity + $scaling_permeabilityinv));
         if [[ $test != 0 ]]; then
   

          echo "Degree                 : $k"
          echo "Mesh family            : $mesh_family"
          echo "Stabilization parameter: $stab_par"
          echo "Test case: solution    : $tcsol"
          echo "Output directory       : $outdir"

          outsubdir=${outdir}/${mesh_family}_k${k}_sv${scaling_viscosity}_sp${scaling_permeabilityinv}
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

          # CREATE DATA FILE FOR LATEX
          echo -e "Deg MeshSize NbCells NbFaces SizeSystem L2ErrorU Rate H1errorU Rate L2errorP Rate EnergyError Rate" > $outsubdir/$errorsfile
          echo -e "TwallVHHOSpace TprocVHHOSpace TwallModel TprocModel TwallSolve TprocSolve" > $outsubdir/$timesfile
          for i in `seq 1 $nbmesh`; 
          do
              Degree=$(awk '/Degree:/ {print $NF}' $outsubdir/results-$i.txt)
              MeshSize=$(awk '/MeshSize:/ {print $NF}' $outsubdir/results-$i.txt)
              NbCells=$(awk '/NbCells:/ {print $NF}' $outsubdir/results-$i.txt)
              NbFaces=$(awk '/NbFaces:/ {print $NF}' $outsubdir/results-$i.txt)
              SizeSystem=$(awk '/SizeSystem:/ {print $NF}' $outsubdir/results-$i.txt)
              L2ErrorU=$(awk '/L2ErrorU:/ {print $NF}' $outsubdir/results-$i.txt)
              H1ErrorU=$(awk '/H1ErrorU:/ {print $NF}' $outsubdir/results-$i.txt)
              L2ErrorP=$(awk '/L2ErrorP:/ {print $NF}' $outsubdir/results-$i.txt)
              EnergyError=$(awk '/EnergyError:/ {print $NF}' $outsubdir/results-$i.txt)
              TwallVHHOSpace=$(awk '/TwallVHHOSpace:/ {print $NF}' $outsubdir/results-$i.txt)
              TprocVHHOSpace=$(awk '/TprocVHHOSpace:/ {print $NF}' $outsubdir/results-$i.txt)
              TwallModel=$(awk '/TwallModel:/ {print $NF}' $outsubdir/results-$i.txt)
              TprocModel=$(awk '/TprocModel:/ {print $NF}' $outsubdir/results-$i.txt)
              TwallSolve=$(awk '/TwallSolve:/ {print $NF}' $outsubdir/results-$i.txt)
              TprocSolve=$(awk '/TprocSolve:/ {print $NF}' $outsubdir/results-$i.txt)
              echo -e "$TwallVHHOSpace $TprocVHHOSpace $TwallModel $TprocModel $TwallSolve $TprocSolve" >> $outsubdir/$timesfile
              if(($i > 1)); then
                imo=$(perl -E "say $i - 1")
                OldMeshSize=$(awk '/MeshSize:/ {print $NF}' $outsubdir/results-$imo.txt)
                OldL2ErrorU=$(awk '/L2ErrorU:/ {print $NF}' $outsubdir/results-$imo.txt)
                L2ErrorURate=$(perl -E "say sprintf(\"%.2f\", log($OldL2ErrorU/$L2ErrorU)/log($OldMeshSize/$MeshSize))")
                OldH1ErrorU=$(awk '/H1ErrorU:/ {print $NF}' $outsubdir/results-$imo.txt)
                H1ErrorURate=$(perl -E "say sprintf(\"%.2f\", log($OldH1ErrorU/$H1ErrorU)/log($OldMeshSize/$MeshSize))")
                OldL2ErrorP=$(awk '/L2ErrorP:/ {print $NF}' $outsubdir/results-$imo.txt)
                L2ErrorPrate=$(perl -E "say sprintf(\"%.2f\", log($OldL2ErrorP/$L2ErrorP)/log($OldMeshSize/$MeshSize))")
                OldEnergyError=$(awk '/EnergyError:/ {print $NF}' $outsubdir/results-$imo.txt)
                EnergyErrorRate=$(perl -E "say sprintf(\"%.2f\", log($OldEnergyError/$EnergyError)/log($OldMeshSize/$MeshSize))")
                echo -e "$Degree $MeshSize $NbCells $NbFaces $SizeSystem $L2ErrorU $L2ErrorURate $H1ErrorU $H1ErrorURate $L2ErrorP $L2ErrorPrate $EnergyError $EnergyErrorRate" >> $outsubdir/$errorsfile
              else
            	  echo -e "$Degree $MeshSize $NbCells $NbFaces $SizeSystem $L2ErrorU -- $H1ErrorU -- $L2ErrorP -- $EnergyError -- " >> $outsubdir/$errorsfile
              fi
          done;
          column -t < $outsubdir/$errorsfile | tee $outsubdir/$errorsfile
          column -t < $outsubdir/$timesfile | tee $outsubdir/$timesfile

          fi # Test pour eviter visc=perm=0

        done # scaling perm

      done # scaling visc

    done # for degree

done # for mesh_family

