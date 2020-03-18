# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Python program to convert .xyz file into mgmol input
#
# use: python coords.xyz > coords.in
#-------------------------------------------------------------------------------
import sys, string

ang2bohr=1.8897269

#read file
ifile=open(sys.argv[1],'r')
lines_of_file=ifile.readlines()

count=0
movable=1 #assuming all atoms can move
dummy=0   #unused flag set to 0

for line in lines_of_file: ## loop over lines of file
  words=string.split(line)
  if len(words)>1:
    if words[0][0:1]!='#':
      name=words[0]+`count`
      x=eval(words[1])*ang2bohr
      y=eval(words[2])*ang2bohr
      z=eval(words[3])*ang2bohr
      
      print name,'\t',dummy,'\t',x,'\t',y,'\t',z,'\t',movable
      count=count+1
