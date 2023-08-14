# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# use: python getBoundingBox.py mgmol_output
#-------------------------------------------------------------------------------
import sys, string
from math import sqrt

input=open(sys.argv[1],'r')

frame=-1
if len(sys.argv)>3:
  frame=eval(sys.argv[3])
  print("Input argument: Frame = {}".format(frame))

lines=input.readlines()

star='*'

##############################################
# count number atoms
def getNumAtoms(lines):
  found_current_line=0
  already_found_one=0
  na=0
  flag=0
  for line in lines: ## loop over lines of file
    if line.count("POSITIONS"):
      flag=1
    if flag==1 and line.count('## '):
      found_current_line=1
      already_found_one =1
      na=na+1
    else:
      found_current_line=0
    if found_current_line!=already_found_one:
      break
  return na
##############################################


na=getNumAtoms(lines)

print("N atoms in file = {}".format(na))


##############################################

searchterm1='POSITIONS'

def getPositions(names,coords,L,fframe):
  na=len(names)
  l=len(L)  ## no of lines in file
  line_min=0  
  cur_frame=0
  for line in range(l): ## loop over lines of file1 
    if line>line_min:
      num_matches1 = L[line].count(searchterm1)
      if num_matches1:
        j=0
        print("Frame: {}, Read positions starting at line {} for {} atoms".format(cur_frame, line, na))
        for line2 in range(line+1,line+na+10):
          if L[line2].count('$$')>0:
            break;
          if L[line2].count('##')>0:
            words=L[line2].split()
            shift=0
            while words[shift]!='##':
              shift=shift+1
            shift=shift+1
            if words[shift]=='*':
              print("found *")
              shift=shift+1
            word=words[shift:]
            name=word[0]
            #print(name)
            x=word[1]
            y=word[2]
            z=word[3]
            names[j]=name
            coords[j]=x+'\t'+y+'\t'+z
            j=j+1
        line_min=line+na+2
        if fframe==cur_frame:
          print("break")
          break
        cur_frame=cur_frame+1

##############################################

coords1=[]
names1=[]
for i in range(0,na):
  coords1.append("")
  names1.append("")

getPositions(names1,coords1,lines,frame)


maxx=-1000.
maxy=-1000.
maxz=-1000.
minx=1000.
miny=1000.
minz=1000.

for i in range(na): 
  word1=coords1[i].split()
  x1=eval(word1[0])
  y1=eval(word1[1])
  z1=eval(word1[2])
  maxx=max(maxx,x1)
  maxy=max(maxy,y1)
  maxz=max(maxz,z1)
  minx=min(minx,x1)
  miny=min(miny,y1)
  minz=min(minz,z1)
  na=na+1

print("N atoms = {}".format(na))
print("Lower left = {}, {}, {}".format(minx,miny,minz))
print("Upper right = {}, {}, {}".format(maxx,maxy,maxz))
