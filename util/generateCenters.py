# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
import sys, string
#run:
#python replicateLRs.py mgmol.cfg nx ny nz

cfg=open(sys.argv[1],'r')

nx=eval(sys.argv[2])
ny=eval(sys.argv[3])
nz=eval(sys.argv[4])

#get lx, ly, lz from cfg file
lx=0.
ly=0.
lz=0.
for line in cfg:
  if "lx" in line:
    word=string.split(line,'=')
    lx=eval(word[1])
  if "ly" in line:
    word=string.split(line,'=')
    ly=eval(word[1])
  if "lz" in line:
    word=string.split(line,'=')
    lz=eval(word[1])

#generate centers with mesh spacing (lx/nx,ly/ny,lz/nz)
hx=lx/nx
hy=ly/ny
hz=lz/nz

#offset
offx=0.5*hx
offy=0.5*hy
offz=0.5*hz

for i in range(nx):
  x=i*hx+offx
  for j in range(ny):
    y=j*hy+offy
    for k in range(nz):
      z=k*hz+offz
      print str(x).ljust(16),str(y).ljust(16),str(z).ljust(16)
