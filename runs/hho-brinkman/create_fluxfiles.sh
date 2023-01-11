#!/bin/bash
#
# Create data file for the convergence of the flux value to the reference flux
#
root=$(pwd)/outputs;
fluxfile=flux_convergence.dat;
foldertemplate=CubeCavityWedge-tets-hexas;
kmax=3;

cd $root

## Grab the reference flux
cd ${foldertemplate}_k$kmax
maxmesh=$(ls results* | wc -w)
refFlux=$(awk '/FluxCavity:/ {print $NF}' results-$maxmesh.txt)
echo -e "\nReference flux: $refFlux (for degree $kmax and mesh $maxmesh)"
echo -e "\n--------------"
cd $root
##

for k in `seq 0 $kmax`; do
  cd ${foldertemplate}_k$k
  echo "Entering ${foldertemplate}_k$k"
  nbmesh=$(ls results* | wc -w)
  
  # Do not count last mesh for max k
  if [ $k -eq $kmax ]; then
    nbmesh=$(($nbmesh-1))
  fi
  
  echo -e "Deg MeshSize NbCells NbFaces SizeSystem RelFluxDiff Rate TwallAssembly TprocAssembly TwallSolve TprocSolve TwallTotal TprocTotal" > $fluxfile
  for i in `seq 1 $nbmesh`; 
  do
      Degree=$(awk '/Degree:/ {print $NF}' results-$i.txt)
      MeshSize=$(awk '/MeshSize:/ {print $NF}' results-$i.txt)
      NbCells=$(awk '/NbCells:/ {print $NF}' results-$i.txt)
      NbFaces=$(awk '/NbFaces:/ {print $NF}' results-$i.txt)
      SizeSystem=$(awk '/SizeSystem:/ {print $NF}' results-$i.txt)
      Flux=$(awk '/FluxCavity:/ {print $NF}' results-$i.txt)
      RelFluxDiff=$(perl -E "say abs($Flux - $refFlux)/abs($refFlux)");
      TwallVHHOSpace=$(awk '/TwallVHHOSpace:/ {print $NF}' results-$i.txt)
      TprocVHHOSpace=$(awk '/TprocVHHOSpace:/ {print $NF}' results-$i.txt)
      TwallModel=$(awk '/TwallModel:/ {print $NF}' results-$i.txt)
      TprocModel=$(awk '/TprocModel:/ {print $NF}' results-$i.txt)
      TwallAssembly=$(perl -E "say $TwallVHHOSpace + $TwallModel");
      TprocAssembly=$(perl -E "say $TprocVHHOSpace + $TprocModel");
      TwallSolve=$(awk '/TwallSolve:/ {print $NF}' results-$i.txt)
      TprocSolve=$(awk '/TprocSolve:/ {print $NF}' results-$i.txt)
      TwallTotal=$(perl -E "say $TwallAssembly + $TwallSolve");
      TprocTotal=$(perl -E "say $TprocAssembly + $TprocSolve");
      if(($i > 1)); then
        imo=$(perl -E "say $i - 1")
        OldMeshSize=$(awk '/MeshSize:/ {print $NF}' results-$imo.txt)
        OldFlux=$(awk '/FluxCavity:/ {print $NF}' results-$imo.txt)
        OldRelFluxDiff=$(perl -E "say abs($OldFlux - $refFlux)/abs($refFlux)");
        FluxRate=$(perl -E "say sprintf(\"%.2f\", log($OldRelFluxDiff/$RelFluxDiff)/log($OldMeshSize/$MeshSize))")
        echo -e "$Degree $MeshSize $NbCells $NbFaces $SizeSystem $RelFluxDiff $FluxRate $TwallAssembly $TprocAssembly $TwallSolve $TprocSolve $TwallTotal $TprocTotal" >> $fluxfile
      else
    	  echo -e "$Degree $MeshSize $NbCells $NbFaces $SizeSystem $RelFluxDiff -- $TwallAssembly $TprocAssembly $TwallSolve $TprocSolve $TwallTotal $TprocTotal" >> $fluxfile
      fi
  done;
  column -t < $fluxfile | tee $fluxfile
  
  echo -e "\n-----------------"
  cd $root
done;


