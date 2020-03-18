# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#usage: python thermalizeInput mgmol.in T > new_mgmol_in
#
import sys, string

from gaussianDistribution import boxmuller
from atomicMasses import getMassAtom
from math import sqrt

#Assign normally distributed velocities

# kb_au=8.617343e-5 [eV/K] / 27.211608 [eV/Ha]
kb_au  = 3.16678939e-06   # [Ha/K]
nmass = 1822.89

# list possible species in a dictionary with their count
species={'H':0,'He':0,'Li':0,'Be':0,'B':0,'C':0,'N':0,'O':0,'F':0,'Na':0,
         'Mg':0,'Al':0,'Si':0,'P':0,'S':0,'Cl':0,'Ca':0,'D':0}
 

input=open(sys.argv[1],'r')
L=input.readlines()

T=eval(sys.argv[2])
kT=kb_au*T

#read file and get atoms
names=[]
for line in L: ## loop over lines of file
  words=string.split(line)
  if len(words)==6 or len(words)==9:
    name1=words[0][0:1]
    name2=words[0][0:2]
    if name2 in species.keys() or name1 in species.keys():
      names.append(words[0]) 

vel=[]
natoms=len(names)
for n in range(0,natoms,2):
  v0,v1=boxmuller()
  
  v2,v3=boxmuller()
  
  v4,v5=boxmuller()
  
  mass=getMassAtom(names[n])*nmass
  #Scale velocites, Maxwell distribution has variance kT/m
  # 0.5*m*<v_i^2>=0.5*kb*T
  std=sqrt(kT/mass)
  vel.append(v0*std)
  vel.append(v1*std)
  vel.append(v2*std)
  
  if( n<natoms-1 ):
    mass=getMassAtom(names[n+1])*nmass
    std=sqrt(kT/mass)

    vel.append(v3*std)
    vel.append(v4*std)
    vel.append(v5*std)


ia=0
for line in L: ## loop over lines of file
  word=string.split(line)
  if len(word)==6 or len(word)==9:
    name1=word[0][0:1]
    name2=word[0][0:2]
    if name2 in species.keys() or name1 in species.keys():
      print word[0].ljust(9),word[1].ljust(3),word[2].ljust(9),\
            word[3].ljust(9),word[4].ljust(9),word[5].ljust(3),
      if eval(word[5])>0:
        print '%9.7f  %9.7f  %9.7f' % (vel[3*ia],vel[3*ia+1],vel[3*ia+2])
      else:
        print '%9.7f  %9.7f  %9.7f' % (0.,0.,0.)
      ia=ia+1
    else:
      print line,
  else:
    print line,
