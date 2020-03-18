# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#source /usr/apps/python/default/setup.csh
#python averageMDpositions.py . > average_positions.dat
import sys, string, os
from numpy import *
from pdb_tools import printPDBfile
from collections import defaultdict

au2ps = 2.418885e-05

times=[]
#coords={}
coords=defaultdict(list)
nsteps=0
time=0.
dt=0.

def appendPositions(filename):
  global coords
  global nsteps
  global time
  global dt
  
  found=0
  found_force=0

  #print '#',filename
  file=open(filename,'r')
  L1=file.readlines()
  for line in L1: ## loop over lines of file
    words=string.split(line)
    if len(words)>1:
      #find positions
      if 'Stepper' in words and 'Forces:' in words:
        found_force=1
        found=0
        nsteps=nsteps+1
        time=time+dt*au2ps
        times.append(time)
      if words[0][0:2]=='##' and found_force:
        name=words[1]
        found=found+1
        valx=eval(words[2])
        valy=eval(words[3])
        valz=eval(words[4])
        if name not in coords:
          coords[name]=list()
        coords[name].append([valx,valy,valz])
    else:
        if found>0:
          found_force=0
                    
  file.close()

#main
filesdir=sys.argv[1]
filenames=os.listdir(filesdir)

inputs=[]
for filename in filenames:
  if 'md_run' in filename:
    inputs.append(filename)

inputs.sort()

#read 'dt'
inputs[0]
file=open(inputs[0],'r')
L1=file.readlines()
for line in L1: ## loop over lines of file
  word=string.split(line)
  if len(word)>1:
    if word[0]=='Timestep':
      dt=eval(word[5])

for filename in inputs:
  appendPositions(filename)

#number of steps/ps
nf=int(ceil(1./(dt*au2ps)))

#sampling length in ps
number_ps=3

xyz=[]
names=[]
movables=[]

#loop over atoms
for c in coords.keys():
  #print coords[c]
  #print len(coords[c])
  npts=min(number_ps*nf,len(coords[c]))

  names.append(c)
  
  #build lists with latest npts values
  x=[]
  y=[]
  z=[]
  for i in range(len(coords[c])-npts,len(coords[c])-1):
    x.append(coords[c][i][0])
    y.append(coords[c][i][1])
    z.append(coords[c][i][2])
  
  #calculate averages with latest npts values
  avex = sum(x[-npts:])/npts
  avey = sum(y[-npts:])/npts
  avez = sum(z[-npts:])/npts

  xyz.append(str(avex)+' '+str(avey)+' '+str(avez))
  movables.append('1')

na=len(names)
printPDBfile(na,names,xyz,movables,'')
