# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#python ~/SVN/MGmol/mgmol/trunk/util/separateCenters.py lrs.in > ! new_lrs.in
import sys, string
import random

#amplitude of random displacement added to degenerate centers
jitter=0.4

random.seed( 11234 )

input_file=open(sys.argv[1],'r')

lines = input_file.readlines()
lines.sort()

def compare(x1,y1,z1,x2,y2,z2):
  d2=1.e-4
  dx=abs(x1-x2)
  dy=abs(y1-y2)
  dz=abs(z1-z2)
  
  if dx*dx+dy*dy+dz*dz<d2:
    return 1
  return 0


######main code############

xold=-100000.
yold=-100000.
zold=-100000.

bx=[]
by=[]
bz=[]
for line in lines:
  words=string.split(line)
  x=eval(words[0])
  y=eval(words[1])
  z=eval(words[2])
  if compare(x,y,z,xold,yold,zold)>0:
    jx=(2.*random.random()-1.0)*jitter
    jy=(2.*random.random()-1.0)*jitter
    jz=(2.*random.random()-1.0)*jitter
    print x+jx,'\t',y+jy,'\t',z+jz
  else:
    print x,'\t',y,'\t',z

  xold=x
  yold=y
  zold=z
