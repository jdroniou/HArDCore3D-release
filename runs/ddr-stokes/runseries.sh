#!/bin/bash
#
# Execute hho file on series of meshes, and calculate outputs
#

# Options
if [[ $1 == "help" ]]; then
    echo -e "\nExecute tests using parameters in data.sh"
    exit;
fi;

# Load data
. data.sh

# File for times
timesfile="times.dat"

# Directories
executable_name="ddr-stokes"
make $executable_name

origin=$(pwd)
if [ ! -f ../directories.sh ]; then
  echo "directories.sh does not exist. Please read the README.txt in the parent folder."
  exit
fi
. ../directories.sh

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
    if($executable -m $meshfile -k $k -s $tcsol --pressure_scaling $pressure_scaling -x $stab_par); then
      # Move outputs
      mv results.txt $outsubdir/results-$i.txt
    fi
done

# CREATE DATA FILE FOR LATEX
echo -e "Deg MeshSize E_HcurlVel Rate1 E_L2GradPre Rate2 E_cHcurlVel Rate3 E_cL2GradPre Rate4 absE_HcurlVel absE_L2GradPre absE_cHcurlVel absE_cL2GradPre" > $outsubdir/$errorsfile
echo -e "TwallDDRCore TprocDDRCore TwallModel TprocModel TwallSolve TprocSolve" > $outsubdir/$timesfile
for i in `seq 1 $nbmesh`; 
do
    Degree=$(awk '/^Degree:/ {print $NF}' $outsubdir/results-$i.txt)
    MeshSize=$(awk '/^MeshSize:/ {print $NF}' $outsubdir/results-$i.txt)
    NbCells=$(awk '/^NbCells:/ {print $NF}' $outsubdir/results-$i.txt)
    NbFaces=$(awk '/^NbFaces:/ {print $NF}' $outsubdir/results-$i.txt)
    NbEdges=$(awk '/^NbEdges:/ {print $NF}' $outsubdir/results-$i.txt)
    DimXCurl=$(awk '/^DimXCurl:/ {print $NF}' $outsubdir/results-$i.txt)
    DimXGrad=$(awk '/^DimXGrad:/ {print $NF}' $outsubdir/results-$i.txt)
    E_HcurlVel=$(awk '/^E_HcurlVel:/ {print $NF}' $outsubdir/results-$i.txt)
    E_L2GradPre=$(awk '/^E_L2GradPre:/ {print $NF}' $outsubdir/results-$i.txt)
    absE_HcurlVel=$(awk '/^absE_HcurlVel:/ {print $NF}' $outsubdir/results-$i.txt)
    absE_L2GradPre=$(awk '/^absE_L2GradPre:/ {print $NF}' $outsubdir/results-$i.txt)
    E_cHcurlVel=$(awk '/^E_cHcurlVel:/ {print $NF}' $outsubdir/results-$i.txt)
    E_cL2GradPre=$(awk '/^E_cL2GradPre:/ {print $NF}' $outsubdir/results-$i.txt)
    absE_cHcurlVel=$(awk '/^absE_cHcurlVel:/ {print $NF}' $outsubdir/results-$i.txt)
    absE_cL2GradPre=$(awk '/^absE_cL2GradPre:/ {print $NF}' $outsubdir/results-$i.txt)
    TwallDDRCore=$(awk '/^TwallDDRCore:/ {print $NF}' $outsubdir/results-$i.txt)
    TprocDDRCore=$(awk '/^TprocDDRCore:/ {print $NF}' $outsubdir/results-$i.txt)
    TwallModel=$(awk '/^TwallModel:/ {print $NF}' $outsubdir/results-$i.txt)
    TprocModel=$(awk '/^TprocModel:/ {print $NF}' $outsubdir/results-$i.txt)
    TwallSolve=$(awk '/^TwallSolve:/ {print $NF}' $outsubdir/results-$i.txt)
    TprocSolve=$(awk '/^TprocSolve:/ {print $NF}' $outsubdir/results-$i.txt)
    echo -e "$TwallDDRCore $TprocDDRCore $TwallModel $TprocModel $TwallSolve $TprocSolve" >> $outsubdir/$timesfile
    if(($i > 1)); then
      imo=$(perl -E "say $i - 1")
      OldMeshSize=$(awk '/^MeshSize:/ {print $NF}' $outsubdir/results-$imo.txt)
      OldE_HcurlVel=$(awk '/^E_HcurlVel:/ {print $NF}' $outsubdir/results-$imo.txt)
      HcurlVelRate=$(perl -E "say sprintf(\"%.2f\", log($OldE_HcurlVel/$E_HcurlVel)/log($OldMeshSize/$MeshSize))")
      OldE_L2GradPre=$(awk '/^E_L2GradPre:/ {print $NF}' $outsubdir/results-$imo.txt)
      L2GradPreRate=$(perl -E "say sprintf(\"%.2f\", log($OldE_L2GradPre/$E_L2GradPre)/log($OldMeshSize/$MeshSize))")
      OldE_cHcurlVel=$(awk '/^E_cHcurlVel:/ {print $NF}' $outsubdir/results-$imo.txt)
      cHcurlVelRate=$(perl -E "say sprintf(\"%.2f\", log($OldE_cHcurlVel/$E_cHcurlVel)/log($OldMeshSize/$MeshSize))")
      OldE_cL2GradPre=$(awk '/^E_cL2GradPre:/ {print $NF}' $outsubdir/results-$imo.txt)
      cL2GradPreRate=$(perl -E "say sprintf(\"%.2f\", log($OldE_cL2GradPre/$E_cL2GradPre)/log($OldMeshSize/$MeshSize))")
      echo -e "$Degree $MeshSize $E_HcurlVel $HcurlVelRate $E_L2GradPre $L2GradPreRate $E_cHcurlVel $cHcurlVelRate $E_cL2GradPre $cL2GradPreRate $absE_HcurlVel $absE_L2GradPre $absE_cHcurlVel $absE_cL2GradPre" >> $outsubdir/$errorsfile
    else
  	  echo -e "$Degree $MeshSize $E_HcurlVel -- $E_L2GradPre -- $E_cHcurlVel -- $E_cL2GradPre -- $absE_HcurlVel $absE_L2GradPre $absE_cHcurlVel $absE_cL2GradPre" >> $outsubdir/$errorsfile
    fi
done;
# For standard output: we do not output all the separate errors for u, p etc.
cat $outsubdir/$errorsfile | cut -d ' ' -f 1-10 | column -t
# Final error files with all errors
column -t < $outsubdir/$errorsfile | tee $outsubdir/$errorsfile > /dev/null
column -t < $outsubdir/$timesfile | tee $outsubdir/$timesfile

