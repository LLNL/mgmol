# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#python plotLRmoves.py md_run0*.out
import sys, string
import matplotlib.pyplot as plt

searchterm1='move'
searchterm2='Orbital'
searchterm3='CONFIGURATION'

#unit conversion factors
au2ps=2.418885e-05

moves=[]
mdsteps=[]
time=0.

for arg in sys.argv[1:]:
  inputfile=open(arg,'r')
  L=inputfile.readlines()
  flag=0

  #read time step
  dt=0
  for line in L:
    num_matches1 = string.count(line, 'Timestep')
    if num_matches1:
      dt=eval(string.split(line)[5])
      break

  #print 'dt=',dt
  #convert to ps
  dt=dt*au2ps

  for line in L:
    num_matches1 = string.count(line, searchterm1)
    num_matches2 = string.count(line, searchterm2)
    if num_matches1 & num_matches2:
      words=string.split(line)
      move=eval(words[4])
      #print move
    num_matches3 = string.count(line, searchterm3)
    num_matches4 = string.count(line, 'WARNING')
    num_matches5 = string.count(line, 'large')
    #get time
    if num_matches3:
      words=string.split(line)
      time=eval(words[1])*dt
    #use intermediate times to plot large moves
    if num_matches4 & num_matches5:
      time=time+0.01*dt
    if num_matches3 or (num_matches4 & num_matches5):
      mdsteps.append(time)
      moves.append(move)
      print time,move

plt.plot(mdsteps,moves,'r.--')
plt.ylabel('Amplitude max. LR move')
plt.xlabel('Time (ps)')

#plt.show()
plt.savefig('delta.png', dpi=100)
