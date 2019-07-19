# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Python program to measure distance between 2 atoms in mgmol output
#
# use: python getDiatomicDistance.py mgmol_output atom1 atom2
#-------------------------------------------------------------------------------
import sys, string
from math import sqrt

input=open(sys.argv[1],'r')
name1=sys.argv[2]
name2=sys.argv[3]

lines=input.readlines()

found1=0
found2=0

#read file and get data for atoms
for line in lines: ## loop over lines of file
  words=line.split()
  if len(words)>0:
    if words[0][0:2]=='##' and words[0][0:3]!='###':
      name=words[1]
      if name[0:1]=='*':
        name=name[1:]
      if( name1==name ):
        x1=eval(words[2])
        y1=eval(words[3])
        z1=eval(words[4])
        found1=1
      if( name2==name ):
        x2=eval(words[2])
        y2=eval(words[3])
        z2=eval(words[4])
        found2=1
      if found1==1 & found2==1:
        d=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
        print ("d[Bohr]={}, d[Ang]={}".format(d, d*0.529177))
        found1=0
        found2=0
