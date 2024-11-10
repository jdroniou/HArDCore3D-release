#!/bin/bash
#
# Execute executable on series of meshes, and calculate outputs
#

# Options
if [[ $1 == "help" ]]; then
    echo -e "\nExecute tests using parameters in data.sh"
    exit;
fi;

# Load data
. data.sh

# File for times and norms
timesfile="times.dat"
normsfile="norms.dat"

# Directories
executable_name="sddr-navier-stokes"
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

outsubdir=${outdir}/${mesh_family}_k${k}_Re${reynolds}_Psca${pressure_scaling}
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
    solutionfile="mesh"$i"_sol"
    # Execute code
    if($executable -m $meshfile -k $k -s $tcsol --boundary-conditions $boundary_conditions -r $reynolds --navier-scaling $navier_scaling --pressure-scaling $pressure_scaling -x $stab_par --solver $solver --newton-verbose $newton_verb --newton-tol $newton_tol --plot $solutionfile); then
      # Move outputs
      mv results.txt $outsubdir/results-$i.txt
      if [ -f $solutionfile"_u.vtu" ]; then
        mv $solutionfile"_u.vtu" $outsubdir
        mv $solutionfile"_p.vtu" $outsubdir
      fi
    fi
done

# CREATE DATA FILE FOR LATEX
echo -e "Deg MeshSize E_HcurlVel Rate1 E_L2Pre Rate2 E_L2GradPre Rate3 E_cHcurlVel Rate4 E_cL2Pre Rate5 E_cL2GradPre Rate6 absE_HcurlVel absE_L2GradPre absE_cHcurlVel absE_cL2GradPre "> $outsubdir/$errorsfile
echo -e "Deg MeshSize discN_HCurlVel discN_HGradPre" > $outsubdir/$normsfile
echo -e "TwallDDRCore TprocDDRCore TwallModel TprocModel TwallSolve TprocSolve" > $outsubdir/$timesfile
for i in `seq 1 $nbmesh`; 
do
    Degree=$(awk '/^Degree:/ {print $NF}' $outsubdir/results-$i.txt)
    MeshSize=$(awk '/^MeshSize:/ {print $NF}' $outsubdir/results-$i.txt)
    NbCells=$(awk '/^NbCells:/ {print $NF}' $outsubdir/results-$i.txt)
    NbFaces=$(awk '/^NbFaces:/ {print $NF}' $outsubdir/results-$i.txt)
    NbEdges=$(awk '/^NbEdges:/ {print $NF}' $outsubdir/results-$i.txt)
    DimSXCurl=$(awk '/^DimSXCurl:/ {print $NF}' $outsubdir/results-$i.txt)
    DimSXGrad=$(awk '/^DimSXGrad:/ {print $NF}' $outsubdir/results-$i.txt)
    E_HcurlVel=$(awk '/^E_HcurlVel:/ {print $NF}' $outsubdir/results-$i.txt)
    E_L2Pre=$(awk '/^E_L2Pre:/ {print $NF}' $outsubdir/results-$i.txt)
    E_L2GradPre=$(awk '/^E_L2GradPre:/ {print $NF}' $outsubdir/results-$i.txt)
    absE_HcurlVel=$(awk '/^absE_HcurlVel:/ {print $NF}' $outsubdir/results-$i.txt)
    absE_L2Pre=$(awk '/^absE_L2Pre:/ {print $NF}' $outsubdir/results-$i.txt)
    absE_L2GradPre=$(awk '/^absE_L2GradPre:/ {print $NF}' $outsubdir/results-$i.txt)
    E_cHcurlVel=$(awk '/^E_cHcurlVel:/ {print $NF}' $outsubdir/results-$i.txt)
    E_cL2Pre=$(awk '/^E_cL2Pre:/ {print $NF}' $outsubdir/results-$i.txt)
    E_cL2GradPre=$(awk '/^E_cL2GradPre:/ {print $NF}' $outsubdir/results-$i.txt)
    absE_cHcurlVel=$(awk '/^absE_cHcurlVel:/ {print $NF}' $outsubdir/results-$i.txt)
    absE_cL2Pre=$(awk '/^absE_cL2Pre:/ {print $NF}' $outsubdir/results-$i.txt)
    absE_cL2GradPre=$(awk '/^absE_cL2GradPre:/ {print $NF}' $outsubdir/results-$i.txt)
    discN_HCurlVel=$(awk '/^discN_HCurlVel:/ {print $NF}' $outsubdir/results-$i.txt)
    discN_HGradPre=$(awk '/^discN_HGradPre:/ {print $NF}' $outsubdir/results-$i.txt)
    TwallDDRCore=$(awk '/^TwallDDRCore:/ {print $NF}' $outsubdir/results-$i.txt)
    TprocDDRCore=$(awk '/^TprocDDRCore:/ {print $NF}' $outsubdir/results-$i.txt)
    TwallModel=$(awk '/^TwallModel:/ {print $NF}' $outsubdir/results-$i.txt)
    TprocModel=$(awk '/^TprocModel:/ {print $NF}' $outsubdir/results-$i.txt)
    TwallSolve=$(awk '/^TwallSolve:/ {print $NF}' $outsubdir/results-$i.txt)
    TprocSolve=$(awk '/^TprocSolve:/ {print $NF}' $outsubdir/results-$i.txt)
    echo -e "$Degree $MeshSize $discN_HCurlVel $discN_HGradPre" >> $outsubdir/$normsfile
    echo -e "$TwallDDRCore $TprocDDRCore $TwallModel $TprocModel $TwallSolve $TprocSolve" >> $outsubdir/$timesfile
    if(($i > 1)); then
      imo=$(perl -E "say $i - 1")
      OldMeshSize=$(awk '/^MeshSize:/ {print $NF}' $outsubdir/results-$imo.txt)
      OldE_HcurlVel=$(awk '/^E_HcurlVel:/ {print $NF}' $outsubdir/results-$imo.txt)
      HcurlVelRate=$(perl -E "say sprintf(\"%.2f\", log($OldE_HcurlVel/$E_HcurlVel)/log($OldMeshSize/$MeshSize))")
      OldE_L2Pre=$(awk '/^E_L2Pre:/ {print $NF}' $outsubdir/results-$imo.txt)
      L2PreRate=$(perl -E "say sprintf(\"%.2f\", log($OldE_L2Pre/$E_L2Pre)/log($OldMeshSize/$MeshSize))")
      OldE_L2GradPre=$(awk '/^E_L2GradPre:/ {print $NF}' $outsubdir/results-$imo.txt)
      L2GradPreRate=$(perl -E "say sprintf(\"%.2f\", log($OldE_L2GradPre/$E_L2GradPre)/log($OldMeshSize/$MeshSize))")
      OldE_cHcurlVel=$(awk '/^E_cHcurlVel:/ {print $NF}' $outsubdir/results-$imo.txt)
      cHcurlVelRate=$(perl -E "say sprintf(\"%.2f\", log($OldE_cHcurlVel/$E_cHcurlVel)/log($OldMeshSize/$MeshSize))")
      OldE_cL2Pre=$(awk '/^E_cL2Pre:/ {print $NF}' $outsubdir/results-$imo.txt)
      cL2PreRate=$(perl -E "say sprintf(\"%.2f\", log($OldE_cL2Pre/$E_cL2Pre)/log($OldMeshSize/$MeshSize))")
      OldE_cL2GradPre=$(awk '/^E_cL2GradPre:/ {print $NF}' $outsubdir/results-$imo.txt)
      cL2GradPreRate=$(perl -E "say sprintf(\"%.2f\", log($OldE_cL2GradPre/$E_cL2GradPre)/log($OldMeshSize/$MeshSize))")
      echo -e "$Degree $MeshSize $E_HcurlVel $HcurlVelRate $E_L2Pre $L2PreRate $E_L2GradPre $L2GradPreRate $E_cHcurlVel $cHcurlVelRate $E_cL2Pre $cL2PreRate $E_cL2GradPre $cL2GradPreRate $absE_HcurlVel $absE_L2GradPre $absE_cHcurlVel $absE_cL2GradPre" >> $outsubdir/$errorsfile
    else
  	  echo -e "$Degree $MeshSize $E_HcurlVel -- $E_L2Pre -- $E_L2GradPre -- $E_cHcurlVel -- $E_cL2Pre -- $E_cL2GradPre -- $absE_HcurlVel $absE_L2GradPre $absE_cHcurlVel $absE_cL2GradPre" >> $outsubdir/$errorsfile
    fi
done;
# For standard output: we do not output all the separate errors for u, p etc.
cat $outsubdir/$errorsfile | cut -d ' ' -f 1-14 | column -t
cat $outsubdir/$normsfile  | column -t
# Final error files with all errors
column -t < $outsubdir/$errorsfile | tee $outsubdir/$errorsfile > /dev/null
column -t < $outsubdir/$normsfile | tee $outsubdir/$normsfile > /dev/null
column -t < $outsubdir/$timesfile | tee $outsubdir/$timesfile

