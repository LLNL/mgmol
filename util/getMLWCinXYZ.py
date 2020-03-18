# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#Create XYZ file with Maximally Localized Wannier Centers
#extracted from MGmol output
import sys, string
from math import sqrt

input=open(sys.argv[1],'r')

lines=input.readlines()

bohr2ang = 0.529177

#count number of centers
searchstring='&&'
n=0
for line in lines:
  num_matches = string.count(line, searchstring)
  if num_matches:
    n=n+1

print n,'\n'

#read centers and convert to Angstroem before writing them out
for line in lines:
  num_matches = string.count(line, searchstring)
  if num_matches:
    words=string.split(line)
    x=eval(words[2])
    y=eval(words[3])
    z=eval(words[4])
    print 'Wa',bohr2ang*x,bohr2ang*y,bohr2ang*z
