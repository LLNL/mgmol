# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#!/usr/bin/python

import sys
import os
from phys import *

if sys.argv[1] == "-h":
    print >>sys.stderr, "Usage: make_neb_images.py coordinate_file1 coordinate_file2 num_images"
    print >>sys.stderr, "       The coordinate files are in mgmol format.

file1=sys.argv[1]
file2=sys.argv[2]
nimages=int(sys.argv[3])

if os.path.isfile(file1):
      f1 = open(file1,"r")
else:
      sys.stderr.write(file1 + " has not been found.\n")
      sys.exit()

if os.path.isfile(file2):
      f2 = open(file2,"r")
else:
      sys.stderr.write(file2 + " has not been found.\n")
      sys.exit()

'''
C1 1 8.57635 9.65453 8.4286 1
C2 1 10.637 7.67104 8.47579 0
H1 2 8.86614 11.1286 9.8714 1
H2 2 8.3191 10.5856 6.58459 1
H3 2 12.4054 8.39866 8.71534 0
H4 2 10.227 6.47019 10.1256 1
H5 2 10.5181 6.43654 6.82989 0
F1 3 6.22384 8.48947 8.99996 1
'''

atoms1=[]
for line in f1.readlines():
  tmp=line.split()
  if len(tmp)!=0:
      atoms1+=[{'symbol':tmp[0],'pseudo':int(tmp[1]),'pos':D3v(float(tmp[2]),float(tmp[3]),float(tmp[4])),'frozen':int(tmp[5])}]
f1.close()

atoms2=[]
for line in f2.readlines():
  tmp=line.split()
  if len(tmp)!=0:
      atoms2+=[{'symbol':tmp[0],'pseudo':int(tmp[1]),'pos':D3v(float(tmp[2]),float(tmp[3]),float(tmp[4])),'frozen':int(tmp[5])}]
f2.close()

natoms1=len(atoms1)
natoms2=len(atoms2)
if natoms1!=natoms2: 
  sys.exit('files have different numbers of atoms')
else:
  natoms=natoms1

for i in range(nimages):
  fout=open('image'+str(i+1),"w")
  for j in range(natoms):
    #if atoms1[j]['frozen']==1:
      print>>fout,atoms1[j]['symbol'],atoms1[j]['pseudo'],atoms1[j]['pos']+(atoms2[j]['pos']-atoms1[j]['pos'])/(nimages+1)*(i+1),atoms1[j]['frozen']
    #else:
    #  print>>fout,atoms1[j]['symbol'],atoms1[j]['pseudo'],atoms1[j]['pos'],atoms1[j]['frozen']
  fout.close()


