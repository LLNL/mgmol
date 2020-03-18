# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#usage
#python ./plotTotalEnergy.py mgmol_output
import sys, string
import matplotlib.pyplot as plt
import numpy as np

#unit conversion factors
au2ps=2.418885e-05
Ha2eV=27.211608

searchterm1='Total'
searchterm2='Energy'
searchterm3='CONFIGURATION'

fig = plt.figure()
ax = fig.add_subplot(111)

colors = ['r', 'b', 'k', 'g', 'c', 'm']
file = 0

print "Opening file."

for arg in sys.argv[1:]:
  inputfile=open(arg,'r')
  L=inputfile.readlines()

  energies=[]
  times=[]

  #read number of atoms
  natoms=0
  for line in L:
    num_matches1 = string.count(line, 'Number')
    num_matches2 = string.count(line, 'ions')
    if num_matches1 & num_matches2:
      natoms=eval(string.split(line)[4])
      break

  #read time step
  dt=0
  for line in L:
    num_matches1 = string.count(line, 'Timestep')
    if num_matches1:
      dt=eval(string.split(line)[5])
      break

  print 'dt=',dt
  #convert to ps
  dt=dt*au2ps
  #read energies and times
  flag=0
  for line in L:
    num_matches3 = string.count(line, searchterm3)
    if num_matches3:
      words=string.split(line)
      time=eval(words[1])*dt
    num_matches1 = string.count(line, searchterm1)
    num_matches2 = string.count(line, searchterm2)
    if num_matches1 & num_matches2:
      if 'Hartree' in line:
          continue
      words=string.split(line)
      energy=eval(words[2])*Ha2eV/natoms
      times.append(time)
      energies.append(energy)
      print 'Energy=',energy

  # plot results
  if file < len(colors):
      color = colors[file]
  else:
      color = 'r'
  ax.plot(times,energies, color+'.--')
  file += 1

  np.savez('energies'+str(file)+'.npz', times, energies)

ax.set_ylabel('Total Energy (eV/atom)')
ax.set_xlabel('time (ps)')

plt.show()
plt.savefig('totalE.png', dpi=100)
