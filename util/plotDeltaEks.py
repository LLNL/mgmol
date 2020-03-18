# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#usage:
#python plotDeltaEks.py mgmol.out

import sys, string
import matplotlib.pyplot as plt

inputfile=open(sys.argv[1],'r')
lines=inputfile.readlines()

deltas=[]

for line in lines:
  num_matches1 = string.count(line, 'delta')
  num_matches2 = string.count(line, '%%')
  if num_matches1 & num_matches2:
    words=string.split(line)
    if len(words)>9:
      deltae=eval(words[9])
      deltas.append(deltae)

plt.plot(deltas,'r.--')
plt.ylabel('delta Eks')
plt.xlabel('inner iterations')
plt.axis([0.,len(deltas),1.e-7,1.e-3])
plt.yscale('log')

plt.show()
#plt.savefig('delta.png', dpi=100)
