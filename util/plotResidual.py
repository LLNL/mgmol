# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
import sys, string
import matplotlib.pyplot as plt

residuals=[]

inputfile=open(sys.argv[1],'r')
L=inputfile.readlines()
flag=0
for line in L:
  num_matches1 = string.count(line, 'Residual')
  num_matches2 = string.count(line, '=')
  if num_matches1 & num_matches2:
    words=string.split(line)
    if len(words)>4:
      residual=eval(words[4])
      residuals.append(residual)
      #print residual

plt.plot(residuals,'r.--')
plt.ylabel('residual')
plt.xlabel('inner iterations')
plt.axis([0.,len(residuals),1.e-3,1.e-1])
plt.yscale('log')

#plt.show()
plt.savefig('residual.png', dpi=100)
