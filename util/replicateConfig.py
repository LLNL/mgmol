# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#run:
#python replicateConfig.py filename.cfg nx ny nz > new_filename.cfg
import sys, string
input=open(sys.argv[1],'r')
nx=eval(sys.argv[2])
ny=eval(sys.argv[3])
nz=eval(sys.argv[4])

Li=input.readlines()
for line in Li: ## loop over lines of file
  word=string.split(line,'=')
  if len(word)>0:
    if word[0]=='nempty':
      print 'nempty=',nx*ny*nz*eval(word[1])
      continue
    if word[0]=='nx':
      print 'nx=',nx*eval(word[1])
      continue
    if word[0]=='ny':
      print 'ny=',ny*eval(word[1])
      continue
    if word[0]=='nz':
      print 'nz=',nz*eval(word[1])
      continue
    if word[0]=='lx':
      print 'lx=',nx*eval(word[1])
      continue
    if word[0]=='ly':
      print 'ly=',ny*eval(word[1])
      continue
    if word[0]=='lz':
      print 'lz=',nz*eval(word[1])
      continue
    if word[0]=='ox':
      print 'ox=',nx*eval(word[1])
      continue
    if word[0]=='oy':
      print 'oy=',ny*eval(word[1])
      continue
    if word[0]=='oz':
      print 'oz=',nz*eval(word[1])
      continue
    if word[0]=='atol':
      print 'atol=',nx*ny*nz*eval(word[1])
      continue
    print line,
