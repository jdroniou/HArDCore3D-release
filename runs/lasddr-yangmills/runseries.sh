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
executable_name="lasddr-yangmills"

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
    if [[ ${#iter[@]} < $nbmesh-1 ]]; then
        j=0
    else
        j=$(($i-1))
    fi
    meshfile=$meshdir"/"$(echo ${mesh[$i]} | cut -d ':' -f 2)
    echo -e "------------------------------------------------------------------------------"
    echo -e "Mesh $i out of $nbmesh: $meshfile"
    echo -e "Output directory: $outsubdir"
    echo -e "------------------------------------------------------------------------------"
    # Execute code
    if($executable -m $meshfile -k $k -s $tcsol -x $stab_par -n ${iter[$j]} -v $stop -f $t_final -t $theta -l $nc -c $compute_cn -i $itersolver -q $constr_IC -b $nl_disc); then
      # Move outputs
      mv results.txt $outsubdir/results-$i.txt
      mv condition-numbers.txt $outsubdir/condition-numbers-$i.txt
    fi
done

# CREATE DATA FILE FOR LATEX
echo -e "Deg MeshSize dt E_L2Elec RateElec E_L2Pot RatePot MaxRes MaxConstr" > $outsubdir/$errorsfile
echo -e "TwallDDRCore TprocDDRCore TwallModel TprocModel TwallAssemble TprocAssemble TwallSolve TprocSolve TwallRun TprocRun" > $outsubdir/$timesfile
for i in `seq 1 $nbmesh`; 
do
    Degree=$(awk '/^Degree:/ {print $NF}' $outsubdir/results-$i.txt)
    MeshSize=$(awk '/^MeshSize:/ {print $NF}' $outsubdir/results-$i.txt)
    Timestep=$(awk '/^Timestep:/ {print $NF}' $outsubdir/results-$i.txt)
    NbCells=$(awk '/^NbCells:/ {print $NF}' $outsubdir/results-$i.txt)
    NbFaces=$(awk '/^NbFaces:/ {print $NF}' $outsubdir/results-$i.txt)
    NbEdges=$(awk '/^NbEdges:/ {print $NF}' $outsubdir/results-$i.txt)
    DimXCurl=$(awk '/^DimXCurl:/ {print $NF}' $outsubdir/results-$i.txt)
    DimXGrad=$(awk '/^DimXGrad:/ {print $NF}' $outsubdir/results-$i.txt)
    E_L2Elec=$(awk '/^E_L2Elec:/ {print $NF}' $outsubdir/results-$i.txt)
    E_L2Pot=$(awk '/^E_L2Pot:/ {print $NF}' $outsubdir/results-$i.txt)
    MaxResidual=$(awk '/^MaxResidual:/ {print $NF}' $outsubdir/results-$i.txt)
    MaxConstraintDiff=$(awk '/^MaxConstraintDiff:/ {print $NF}' $outsubdir/results-$i.txt)
    TwallDDRCore=$(awk '/^TwallDDRCore:/ {print $NF}' $outsubdir/results-$i.txt)
    TprocDDRCore=$(awk '/^TprocDDRCore:/ {print $NF}' $outsubdir/results-$i.txt)
    TwallModel=$(awk '/^TwallModel:/ {print $NF}' $outsubdir/results-$i.txt)
    TprocModel=$(awk '/^TprocModel:/ {print $NF}' $outsubdir/results-$i.txt)
    TwallAssemble=$(awk '/^TwallAssemble:/ {print $NF}' $outsubdir/results-$i.txt)
    TprocAssemble=$(awk '/^TprocAssemble:/ {print $NF}' $outsubdir/results-$i.txt)
    TwallSolve=$(awk '/^TwallSolve:/ {print $NF}' $outsubdir/results-$i.txt)
    TprocSolve=$(awk '/^TprocSolve:/ {print $NF}' $outsubdir/results-$i.txt)
    TwallRun=$(awk '/^TwallRun:/ {print $NF}' $outsubdir/results-$i.txt)
    TprocRun=$(awk '/^TprocRun:/ {print $NF}' $outsubdir/results-$i.txt)
    echo -e "$TwallDDRCore $TprocDDRCore $TwallModel $TprocModel $TwallAssemble $TprocAssemble $TwallSolve $TprocSolve $TwallRun $TprocRun" >> $outsubdir/$timesfile
    if(($i > 1)); then
      imo=$(perl -E "say $i - 1")
      OldMeshSize=$(awk '/^MeshSize:/ {print $NF}' $outsubdir/results-$imo.txt)
      OldE_L2Elec=$(awk '/^E_L2Elec:/ {print $NF}' $outsubdir/results-$imo.txt)
      OldE_L2Pot=$(awk '/^E_L2Pot:/ {print $NF}' $outsubdir/results-$imo.txt)
      E_L2ElecRate=$(perl -E "say sprintf(\"%.2f\", log($OldE_L2Elec/$E_L2Elec)/log($OldMeshSize/$MeshSize))")
      E_L2PotRate=$(perl -E "say sprintf(\"%.2f\", log($OldE_L2Pot/$E_L2Pot)/log($OldMeshSize/$MeshSize))")
      echo -e "$Degree $MeshSize $Timestep $E_L2Elec $E_L2ElecRate $E_L2Pot $E_L2PotRate $MaxResidual $MaxConstraintDiff" >> $outsubdir/$errorsfile
    else
  	  echo -e "$Degree $MeshSize $Timestep $E_L2Elec -- $E_L2Pot -- $MaxResidual $MaxConstraintDiff" >> $outsubdir/$errorsfile
    fi
done;
# For standard output: we do not output all the separate errors for u, p etc.
cat $outsubdir/$errorsfile | cut -d ' ' -f 1-10 | column -t
# Final error files with all errors
column -t < $outsubdir/$errorsfile | tee $outsubdir/$errorsfile > /dev/null
column -t < $outsubdir/$timesfile | tee $outsubdir/$timesfile

