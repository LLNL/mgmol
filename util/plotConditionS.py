# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
import sys, string
import matplotlib.pyplot as plt

searchterm1='Condition'
searchterm2='CONFIGURATION'

conds=[]
mdsteps=[]

inputfile=open(sys.argv[1],'r')
L=inputfile.readlines()
flag=0
for line in L:
  num_matches1 = string.count(line, searchterm1)
  if num_matches1:
    words=string.split(line)
    cond=eval(words[4])
  num_matches2 = string.count(line, searchterm2)
  if num_matches2:
    words=string.split(line)
    it=eval(words[1])
    mdsteps.append(it)
    conds.append(cond)

plt.plot(mdsteps,conds,'r.--')
plt.ylabel('condition S')
plt.xlabel('MD step')

plt.show()
#plt.savefig('delta.png', dpi=100)
