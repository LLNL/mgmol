# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#run:
#python replicateLRs.py coords.in mgmol_quench.cfg 2 2 2
import sys, string
lrs=open(sys.argv[1],'r')
cfg=open(sys.argv[2],'r')

nx=eval(sys.argv[3])
ny=eval(sys.argv[4])
nz=eval(sys.argv[5])

#get lx, ly, lz from config file
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

# loop over lines of LRs file and print replicated centers
for line in lrs:
  word=string.split(line)
  if len(word)>0:
    for i in range(nx):
      x=eval(word[0])+i*lx
      for j in range(ny):
        y=eval(word[1])+j*ly
        for k in range(nz):
          z=eval(word[2])+k*lz
          print str(x).ljust(16),str(y).ljust(16),str(z).ljust(16)
