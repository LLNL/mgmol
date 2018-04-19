#!/bin/csh
#
#Remove sections surrounded by ifdef/endif in src directory
#and create a new directory with modified files
#
#Requires executable 'unifdef'
#
set dirs= ( src src/tools src/pb src/DistMatrix src/sparse_linear_algebra src/linear_algebra src/numerical_kernels src/local_matrices src/radial)
foreach dir ($dirs)
  echo $dir
  set newdir = new_$dir
  if(! -d $newdir)mkdir $newdir
  cp $dir/module.mk $newdir
  foreach x ($dir/*.cc $dir/*.h)
     echo $x, $newdir
     ~/bin/unifdef -U__MGMOL_DAVIDSON__ $x > new_$x
  end
end
cp src/Makefile new_src


