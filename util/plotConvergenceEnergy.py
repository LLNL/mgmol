# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#usage: python plotConvergenceEnergy.py mgmol_output
import sys, string
import matplotlib.pyplot as plt

conv_energy=10000.

markers=['r.--','b.--','g.--']

i=0
for filename in sys.argv[1:]:
  energies=[]

  inputfile=open(filename,'r')
  lines=inputfile.readlines()

  flag=0
  na=0
  for line in lines:
    if line.count('Number of ions'):
      words=line.split()
      na=eval(words[4])
      print('na = {}'.format(na))
    if line.count('ENERGY') & line.count('%%'):
      words=line.split()
      energy=eval(words[5][:-1])
      energies.append(energy)
      if conv_energy>energy:
        conv_energy=energy

  print('Reference energy [Ha/atom] = {}'.format(conv_energy/na))
  deltaes=[]
  for energy in energies:
    deltaes.append((energy-conv_energy)/na)

  plt.plot(deltaes,markers[i])
  plt.axis([0.,len(deltaes),10.*deltaes[-2],deltaes[0]])
  i=i+1

plt.ylabel('error Eks/atom [Ha]', fontsize=12)
plt.xlabel('outer iterations', fontsize=12)
plt.yscale('log')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

#plt.show()
plt.savefig('errorEnergy.png', dpi=100)
