# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#plot trajectory of an atom from data in file
#usage:
#  python plotTrajectory.py filename atom_name

import sys, string
import matplotlib.pyplot as plt

inputfile=open(sys.argv[1],'r')
aname    =sys.argv[2]

xx=[]
yy=[]
zz=[]

lines=inputfile.readlines()
flag=0
for line in lines:
  num_matches = string.count(line, aname)
  if num_matches:
    words=string.split(line)
    #print words
    x=eval(words[2])
    y=eval(words[3])
    z=eval(words[4])
    xx.append(x)
    yy.append(y)
    zz.append(z)

if len(xx)==0:
  print '\nERROR: Atom '+aname+' not found in file '+sys.argv[1]+'\n'
  sys.exit(1)
  
#subtract averages
avg=0.
for x in xx:
  avg=avg+x
avg=avg/len(xx)
for i in range(len(xx)):
  xx[i]=xx[i]-avg
  
avg=0.
for y in yy:
  avg=avg+y
avg=avg/len(yy)
for i in range(len(yy)):
  yy[i]=yy[i]-avg
  
avg=0.
for z in zz:
  avg=avg+z
avg=avg/len(zz)
for i in range(len(zz)):
  zz[i]=zz[i]-avg

#plot results
plt.plot(xx,'r.--')
plt.plot(yy,'b.--')
plt.plot(zz,'c.--')
plt.ylabel('x, y, z')
plt.xlabel('MD step')

plt.show()
#plt.savefig('delta.png', dpi=100)
