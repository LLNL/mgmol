#!/bin/csh
echo '#Energies:'
foreach file ( `ls -lt Cl2*.out| awk '{print $9}'` )
  set energy = `grep %% $file |tail -n1|awk '{split ($6,a,","); print a[1]}'`
  set d = `grep Cl2 $file | grep '##' | awk '{ print $5}'`
  echo $d  $energy
end
