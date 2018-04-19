#!/bin/csh
set dirnumber = `pwd | awk '{split($1,a,"_");print a[2]}' | awk '{split($1,a,"/");print a[1]}'`
set dirname=ACHE/CUT5/D1/D1_$dirnumber

#purge quench files
foreach file ( `ls -lt quench*_1.out| awk '{print $9}'` )
  echo $file ...
  set dir=`grep snapshot $file | tail -n6 | grep close| awk '{split($4,a,"/");print a[1]}'`
  echo 'remove ' $dir ...
  rm -rf $dir
end

foreach file ( `ls -lt quench*_2.out| awk '{print $9}'` )
  echo $file ...
  set dir=`grep snapshot $file | tail -n6 | grep close| awk '{split($4,a,"/");print a[1]}'`
  echo 'remove ' $dir ...
  rm -rf $dir
end

#purge intermediate files
foreach file ( `ls -lt md*.out| awk '{print $9}'` )
  echo $file ...
  set dir=`grep snapshot $file | grep close| awk '{split($4,a,"/");print a[1]}'`
  echo $dir ...
  if ( $#dir>1 ) then
    echo 'remove ' $dir[1]
    rm -rf $dir[1]
  endif
end

#purge old files
foreach file ( `ls -lt md*.out| awk '{print $9}'` )
  echo $file ...
  set no = `echo $file | awk '{split($1,a,"d");print a[2]}' | awk '{split($1,a,".");print a[1]}'`
  set r = `expr $no % 10`
  if( ($r>0) && ($r != 5) )then
    set dir=`grep snapshot $file | tail -n6 | grep close| awk '{split($4,a,"/");print a[1]}'`
    set dirno1 = `echo $dir | awk '{split($1,a,"_");print a[2]}'`
    set dirno2 = `echo $dir | awk '{split($1,a,"_");print a[3]}'`
    if( $dirno1<11 )then
      if( $dirno2<370 )then
        echo 'remove ' $dir
        rm -rf $dir
      endif
    endif
  endif
end
