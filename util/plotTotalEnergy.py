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

fig = plt.figure()
ax = fig.add_subplot(111)

colors = ['r', 'b', 'k', 'g', 'c', 'm']
ifile = 0

print("Opening file...")

for arg in sys.argv[1:]:
  inputfile=open(arg,'r')
  lines = inputfile.readlines()

  totalE=[]
  kineticE=[]
  times=[]

  #read number of atoms
  natoms = 0
  for line in lines:
    if line.count('Number') & line.count('ions'):
      natoms=eval(line.split()[4])
      print("natoms = {}".format(natoms))
      break

  #read time step
  dt=0
  for line in lines:
    if line.count('Timestep'):
      dt=eval(line.split()[5])
      break

  print("dt = {}".format(dt))
  #convert to ps
  dt=dt*au2ps

  #read totalE and times
  flag=0
  for line in lines:
    if line.count('CONFIGURATION'):
      words=line.split()
      time=eval(words[1])*dt
    if line.count('Total') & line.count('Energy'):
      if 'Hartree' in line:
          continue
      words=line.split()
      energy=eval(words[2])*Ha2eV/natoms
      times.append(time)
      totalE.append(energy)
      print("Energy = {}".format(energy))

  # plot results
  print("Plot results...")
  print("Number of energy values: {}".format(len(totalE)))
  if ifile < len(colors):
      color = colors[ifile]
  else:
      color = 'r'
  ax.plot(times,totalE, color+'.--')
  ifile += 1

  np.savez('totalE'+str(ifile)+'.npz', times, totalE)

ax.set_ylabel('Energy (eV/atom)')
ax.set_xlabel('time (ps)')

#turn on either option to show or save into a file
#plt.show()
plt.savefig('totalE.png',bbox_inches='tight')
