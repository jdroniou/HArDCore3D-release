#!/bin/bash
#
# Execute hho-brinkman on series of meshes, and calculate outputs
#

# If we pass an argument, it's either "help", or anything else in which case it means we run the script just to
# generate the data file (we do not re-run the executables)
runexec="true";
if (( $# > 0 )); then
  if [[ $1 == "help" ]]; then
    echo -e "\nExecute tests using parameters in data.sh. If any parameter (other than 'help') is provided, do not re-run the tests but only create the data file."
    exit;
  else
    runexec="false";
    echo -e "\nJust creating data file using existing result files...\n"
  fi;
fi;

# Load data
. data.sh

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

outsubdir=${outdir}/${mesh_family}_k${k}_sv${scaling_viscosity}_sp${scaling_permeabilityinv}

echo "Degree                 : $k"
echo "Mesh family            : $mesh_family"
if [ $runexec == "true" ]; then
  echo "Stabilization parameter: $stab_par_stokes / $stab_par_darcy"
  echo "Test case: solution    : $tcsol"
fi
echo -e "Output directory       : $outsubdir\n"

###
# EXECUTE FOR EACH MESH
# Do not touch the results directory if we do not rerun the exec
if [ $runexec == "true" ]; then

  if [ ! -d $outsubdir ]; then
      mkdir -p $outsubdir
  else
      \rm -r $outsubdir/*
  fi

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

fi # $runexec == true

# CREATE DATA FILE FOR LATEX
echo -e "Deg MeshSize NbCells NbFaces SizeSystem L2ErrorU Rate H1errorU Rate L2errorP Rate EnergyError Rate" > $outsubdir/$errorsfile
echo -e "TwallVHHOSpace TprocVHHOSpace TwallModel TprocModel TwallSolve TprocSolve" > $outsubdir/$timesfile

# Check how many result files are available
nbresults=$(ls $outsubdir/results* | wc -w)

for i in `seq 1 $nbresults`; 
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
echo -e ""
column -t < $outsubdir/$timesfile | tee $outsubdir/$timesfile

