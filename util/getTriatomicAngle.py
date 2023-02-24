# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Python program to measure angle between 3 atoms in mgmol output
#
# use: python getTriatomicAngle.py mgmol_output atom1 atom2 atom3
#-------------------------------------------------------------------------------
import sys, string
from math import sqrt, pi, acos

input=open(sys.argv[1],'r')
name1=sys.argv[2]
name2=sys.argv[3]
name3=sys.argv[4]

lines=input.readlines()

found1=0
found2=0
found3=0

#read file and get data for atoms
for line in lines: ## loop over lines of file
  words=line.split()
  if len(words)>1:
    if words[0][0:2]=='##':
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
      if( name3==name ):
        x3=eval(words[2])
        y3=eval(words[3])
        z3=eval(words[4])
        found3=1
      if found1==1 & found2==1 & found3==1:
        d12=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
        d13=sqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
        d23=sqrt((x2-x3)**2+(y2-y3)**2+(z2-z3)**2)
        print("d12 [Bohr]={}, d12 [Ang]={}".format(d12, d12*0.529177))
        print("d13 [Bohr]={}, d13 [Ang]={}".format(d13, d13*0.529177))
        print("d23 [Bohr]={}, d23 [Ang]={}".format(d23, d23*0.529177))
        
        v12_x=(x2-x1)/d12
        v12_y=(y2-y1)/d12
        v12_z=(z2-z1)/d12
        
        v13_x=(x3-x1)/d13
        v13_y=(y3-y1)/d13
        v13_z=(z3-z1)/d13
        
        angle = acos( v12_x*v13_x+v12_y*v13_y+v12_z*v13_z )
        print("angle={}".format(angle*180./pi))
        
        found1=0
        found2=0
        found3=0
