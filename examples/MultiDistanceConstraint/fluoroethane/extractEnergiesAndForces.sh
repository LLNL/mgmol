#!/bin/csh
set utildir = ~/SVN/MGmol/mgmol/trunk/util

set atom1=C1
set atom2=F3
set atom3=C2
set atom4=H7
set atom5=F3
set atom6=H7
set alpha1=1.
set alpha2=1.
set alpha3=-1.

echo '# distance  energy  force'

foreach filename (`ls -l *.out| awk '{print $9}'`)
  #echo $filename
  set force=`grep 'Max. |' $filename |tail -n1|awk '{print $4}'`
  #echo $force
  set energy=`grep 'IONIC CONF' $filename |tail -n1|awk '{print $7}'|cut -c1-11`
  set nc=`grep constraint $filename|grep 'force =' |wc|awk '{print $1}'`
  set n = `grep constraint $filename|grep 'force ='|wc -l`
  if( $n>0 )then
    set forcec = `grep constraint $filename|grep 'force =' |tail -n1|awk '{print $15}'`
  else
    set forcec = 0.0
  endif
  set string1=`python $utildir/getDiatomicDistance.py $filename $atom1 $atom2 |tail -n1|awk '{print $3}'|cut -c1-5`
  set string2=`python $utildir/getDiatomicDistance.py $filename $atom3 $atom4 |tail -n1|awk '{print $3}'|cut -c1-5`
  set string3=`python $utildir/getDiatomicDistance.py $filename $atom5 $atom6 |tail -n1|awk '{print $3}'|cut -c1-5`
  set d=`echo $alpha1 \* $string1+$alpha2 \* $string2+$alpha3 \* $string3|bc -l`
  echo $d"   " $energy"  " $forcec"   #Max. |F|="$force
end
