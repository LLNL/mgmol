# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#how to run in local directory:
#python averageDiffDistances.py . 1. H3041 O8203 -1. H3041 N6678 > distances.dat

import sys, string, os
from math import sqrt, acos
from numpy import *
import matplotlib.pyplot as plt


kb_au=3.16678939e-06 # [Ha/K]
au2ps = 2.418885e-05


distances=[]
times=[]
time=0.


def getMassAtom(name):
  one=name[0]
  mass=0.
  if one[0:2] in spmass.keys():
     mass=spmass[one[0:0+2]]
  else:
     mass=spmass[one[0]]
  return mass

def analyzeDistance(filename,name0,name1,name2,name3,alpha0,alpha1,dt):
  global distances
  global times
  global time
  
  found0=0
  found1=0
  found2=0
  found3=0

  file=open(filename,'r')
  L1=file.readlines()
  for line in L1: ## loop over lines of file
    words=string.split(line)
    if len(words)>1:
      if words[0][0:2]=='##' and words[0][0:3]!='###':
        name=words[1]
        if name[0:1]=='*':
          name=name[1:]

        if( name0==name ):
          x0=eval(words[2])
          y0=eval(words[3])
          z0=eval(words[4])
          found0=1
        if( name1==name ):
          x1=eval(words[2])
          y1=eval(words[3])
          z1=eval(words[4])
          found1=1
        if( name2==name ):
          x2=eval(words[2])
          y2=eval(words[3])
          z2=eval(words[4])
          found2=1
        if( name3==name ):
          x3=eval(words[2])
          y3=eval(words[3])
          z3=eval(words[4])
          found3=1
        if found0==1 & found1==1 & found2==1 & found3==1 :
          d01=sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
          d23=sqrt((x3-x2)**2+(y3-y2)**2+(z3-z2)**2)
          d=eval(alpha0)*d01+eval(alpha1)*d23
          #print "d [Bohr]=",d, ", d [Ang]=",d*0.529177

          times.append(time)
          distances.append(d)
      
          time=time+dt*au2ps

          found0=0
          found1=0
          found2=0
          found3=0
      
  file.close()

#main
filesdir=sys.argv[1]
filenames=os.listdir(filesdir)

alpha0=sys.argv[2]
name0 =sys.argv[3]
name1 =sys.argv[4]

alpha1=sys.argv[5]
name2 =sys.argv[6]
name3 =sys.argv[7]

inputs=[]
for filename in filenames:
  if 'md_run' in filename:
    inputs.append(filename)

inputs.sort()

inputs[0]
file=open(inputs[0],'r')
L1=file.readlines()
for line in L1: ## loop over lines of file
  word=string.split(line)
  if len(word)>1:
    if word[0]=='Timestep':
      dt=eval(word[5])

for filename in inputs:
  analyzeDistance(filename,name0,name1,name2,name3,alpha0,alpha1,dt)

nf=int(ceil(1./(dt*au2ps)))

skip=5
distancep = [ distances[i] for i in range(0, len(distances), skip)]
timesp    = [ times[i]     for i in range(0, len(distances), skip)]

skip=25
aves1 = [sum(distances[i-nf:i])/nf       for i in range(nf,   len(distances), skip)]
aves2 = [sum(distances[i-2*nf:i])/(2*nf) for i in range(2*nf, len(distances), skip)]
aves3 = [sum(distances[i-3*nf:i])/(3*nf) for i in range(3*nf, len(distances), skip)]

times1 = [ times[i] for i in range(nf, len(distances), skip)]
times2 = [ times[i] for i in range(2*nf, len(distances), skip)]
times3 = [ times[i] for i in range(3*nf, len(distances), skip)]

print '#time    distance'
for i in range(len(distances)):
  print times[i],distances[i]

if len(aves1)>1:
  print '#Running average over 1 ps, last value: ',aves1[-1]
if len(aves2)>1:
  print '#Running average over 2 ps, last value: ',aves2[-1]
if len(aves3)>1:
  print '#Running average over 3 ps, last value: ',aves3[-1]


xmax=len(distances)*dt*au2ps
ymin=min(distances)
ymax=max(distances)

plt.figure(1)
plt.subplot(211)
plt.axis([0.,xmax,ymin,ymax])
plt.plot(timesp, distancep, 'go')
plt.ylabel('distance (Bohr)')
#plt.show()

#plt.subplot(212)
#ymin=min(aves2)
#ymax=max(aves2)
#plt.axis([0.,xmax,ymin,ymax])
#plt.plot(times1, aves1, 'ro')
plt.plot(times2, aves2, 'ro')
plt.plot(times3, aves3, 'bo')
#plt.show()
plt.savefig('aves.png')

