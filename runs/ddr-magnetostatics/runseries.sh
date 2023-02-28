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
executable_name="ddr-magnetostatics"
#make $executable_name

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
    if($executable -m $meshfile -k $k -s $tcsol -x $stab_par); then
      # Move outputs
      mv results.txt $outsubdir/results-$i.txt
    fi
done

# CREATE DATA FILE FOR LATEX
echo -e "Deg MeshSize NbCells NbFaces NbEdges DimXCurl DimXDiv SizeSCsystem HcurlHdivError Rate" > $outsubdir/$errorsfile
echo -e "TwallDDRCore TprocDDRCore TwallModel TprocModel TwallSolve TprocSolve" > $outsubdir/$timesfile
for i in `seq 1 $nbmesh`; 
do
    Degree=$(awk '/Degree:/ {print $NF}' $outsubdir/results-$i.txt)
    MeshSize=$(awk '/MeshSize:/ {print $NF}' $outsubdir/results-$i.txt)
    NbCells=$(awk '/NbCells:/ {print $NF}' $outsubdir/results-$i.txt)
    NbFaces=$(awk '/NbFaces:/ {print $NF}' $outsubdir/results-$i.txt)
    NbEdges=$(awk '/NbEdges:/ {print $NF}' $outsubdir/results-$i.txt)
    DimXCurl=$(awk '/DimXCurl:/ {print $NF}' $outsubdir/results-$i.txt)
    DimXDiv=$(awk '/DimXDiv:/ {print $NF}' $outsubdir/results-$i.txt)
    SizeSCsystem=$(awk '/SizeSCsystem:/ {print $NF}' $outsubdir/results-$i.txt)
    EnergyError=$(awk '/EnergyError:/ {print $NF}' $outsubdir/results-$i.txt)
    HcurlHdivError=$(awk '/HcurlHdivError:/ {print $NF}' $outsubdir/results-$i.txt)
    TwallDDRCore=$(awk '/TwallDDRCore:/ {print $NF}' $outsubdir/results-$i.txt)
    TprocDDRCore=$(awk '/TprocDDRCore:/ {print $NF}' $outsubdir/results-$i.txt)
    TwallModel=$(awk '/TwallModel:/ {print $NF}' $outsubdir/results-$i.txt)
    TprocModel=$(awk '/TprocModel:/ {print $NF}' $outsubdir/results-$i.txt)
    TwallSolve=$(awk '/TwallSolve:/ {print $NF}' $outsubdir/results-$i.txt)
    TprocSolve=$(awk '/TprocSolve:/ {print $NF}' $outsubdir/results-$i.txt)
    echo -e "$TwallDDRCore $TprocDDRCore $TwallModel $TprocModel $TwallSolve $TprocSolve" >> $outsubdir/$timesfile
    if(($i > 1)); then
      imo=$(perl -E "say $i - 1")
      OldMeshSize=$(awk '/MeshSize:/ {print $NF}' $outsubdir/results-$imo.txt)
      OldEnergyError=$(awk '/EnergyError:/ {print $NF}' $outsubdir/results-$imo.txt)
      EnergyErrorRate=$(perl -E "say sprintf(\"%.2f\", log($OldEnergyError/$EnergyError)/log($OldMeshSize/$MeshSize))")
      OldHcurlHdivError=$(awk '/HcurlHdivError:/ {print $NF}' $outsubdir/results-$imo.txt)
      HcurlHdivErrorRate=$(perl -E "say sprintf(\"%.2f\", log($OldHcurlHdivError/$HcurlHdivError)/log($OldMeshSize/$MeshSize))")
      echo -e "$Degree $MeshSize $NbCells $NbFaces $NbEdges $DimXCurl $DimXDiv $SizeSCsystem $HcurlHdivError $HcurlHdivErrorRate" >> $outsubdir/$errorsfile
    else
  	  echo -e "$Degree $MeshSize $NbCells $NbFaces $NbEdges $DimXCurl $DimXDiv $SizeSCsystem $HcurlHdivError -- " >> $outsubdir/$errorsfile
    fi
done;
column -t < $outsubdir/$errorsfile | tee $outsubdir/$errorsfile
column -t < $outsubdir/$timesfile | tee $outsubdir/$timesfile

