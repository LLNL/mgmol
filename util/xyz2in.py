# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Python program to convert .xyz file into mgmol input coordinates file
# Optional arguments: [lx ly lz] to define box size and map all atoms into
#                     box (0,0,0)-(lx,ly,lz) using periodic boundary conditions
#
# use: python coords.xyz [lx ly lz] > coords.in
#-------------------------------------------------------------------------------
import sys, string

ang2bohr=1.8897269

#read file
ifile=open(sys.argv[1],'r')

lx = 0.
ly = 0.
lz = 0.
if( len(sys.argv) > 2 ):
  lx = eval(sys.argv[2])
  ly = eval(sys.argv[3])
  lz = eval(sys.argv[4])


lines=ifile.readlines()

count=0
movable=1 #assuming all atoms can move
dummy=0   #unused flag set to 0

for line in lines: ## loop over lines of file
  if count>1:
    words=line.split()
    if len(words)>1:
      name=words[0]+str(count-2)
      x=eval(words[1])*ang2bohr
      y=eval(words[2])*ang2bohr
      z=eval(words[3])*ang2bohr

      if lx > 0.:
        if x<0:
          x = x +lx
        if x>lx:
          x = x -lx

      if ly > 0.:
        if y<0:
          y = y +ly
        if y>ly:
          y = y -ly

      if lz > 0.:
        if z<0:
          z = z +lz
        if z>lz:
          z = z -lz

      print(name,'\t',dummy,'\t',x,'\t',y,'\t',z,'\t',movable)
  count=count+1
