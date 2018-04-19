#!/bin/csh
echo '#Forces:'
foreach file ( `ls -lt Cl2*.out| awk '{print $9}'` )
  grep Cl2 $file | grep '##' | awk '{ print $5, $8 }'
end
