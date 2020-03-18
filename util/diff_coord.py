# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Python program to compare coordinates between 2 mgmol outputs
#
# use: python diff_coord.py mgmol_output1 mgmol_output2
#-------------------------------------------------------------------------------
import sys, string
from math import sqrt

input1=open(sys.argv[1],'r')
input2=open(sys.argv[2],'r')

L1=input1.readlines()
l1=len(L1)  ## no of lines in file
L2=input2.readlines()
l2=len(L2)  ## no of lines in file



##############################################
# count number atoms
def getNumAtoms(L):
  searchterm='## '
  found_current_line=0
  already_found_one=0
  na=0
  for line in L: ## loop over lines of file 
    num_matches = string.count(line, searchterm)
    if num_matches==1:
      found_current_line=1
      already_found_one =1
      na=na+1
    else:
      found_current_line=0
    if found_current_line!=already_found_one:
      break
  return na
##############################################


na1=getNumAtoms(L1)
na2=getNumAtoms(L2)

print 'N atoms in file1=', na1
print 'N atoms in file2=', na2


##############################################
def getCoords(names,coords,L):
  searchterm='## '
  na=len(names)
  l=len(L)  ## no of lines in file
  line_min=0  
  for line in range(l): ## loop over lines of file
    if line>line_min:
      num_matches = string.count(L[line], searchterm)
      if num_matches==1:
        j=0
        print 'Read coords starting at line ', line
        for line2 in range(line,line+na):
          words=string.split(L[line2])
          shift=0
          while words[shift]!='##':
            shift=shift+1
          shift=shift+1
          word=words[shift:]
          name=word[0]
          if name=='*':
            name=word[1]
            word=word[1:]
          if name[0]=='*':
            name=name[1:]
          if name[0]=='D':
            name='H'+name[1:]
          x=word[1]
          y=word[2]
          z=word[3]
          names[j]=name
          coords[j]=x+'\t'+y+'\t'+z
          j=j+1
        line_min=line+na+1
##############################################

coords1=[]
names1=[]
for i in range(0,na1):
  coords1.append(0)
  names1.append(0)

coords2=[]
names2=[]
for i in range(0,na2):
  coords2.append(0)
  names2.append(0)

getCoords(names1,coords1,L1)
getCoords(names2,coords2,L2)

mindr=100.
maxdr=0.
avg=0.
avgx=0.
avgy=0.
avgz=0.
drr=[]
bin=[]
for i in range(0,10):
  bin.append(0)

list_distances=[]
na=0
nfrozen=0
for i in range(na1): 
  word1=string.split(coords1[i])
  for j in range(na2): 
    word2=string.split(coords2[j])
    if names1[i]==names2[j]:
      dx=eval(word1[0])-eval(word2[0])
      dy=eval(word1[1])-eval(word2[1])
      dz=eval(word1[2])-eval(word2[2])
      dr=sqrt(dx*dx+dy*dy+dz*dz)
      drr.append(dr)
      
      list_distances.append(`dr`+' '+names1[i])
      avg=avg+dr
      avgx=avgx+dx
      avgy=avgy+dy
      avgz=avgz+dz
      if dr>maxdr:
        maxdr=dr
        imax=i
        jmax=j
      if dr<mindr:
        mindr=dr
      if dr<1.e-8:
        nfrozen=nfrozen+1
      na=na+1
      print names1[i],'\t',dr

ang2bohr=1.8897269
bohr2ang=1./ ang2bohr

list_distances.sort()
print '-------------------------------------'
print 'List atoms with largest displacements:'
print 'Number of comparisons: ', len(list_distances)
for j in range(min(30,na1)):
  word=string.split(list_distances[na1-j-1])
  print word[1],eval(word[0])*bohr2ang,'[Ang]'

print '-------------------------------------'
print 'Avg. d=',avg/(3*na)*bohr2ang,'[Ang]'
print 'Max. d=',maxdr*bohr2ang,'[Ang]'
print 'Nb. frozen atoms:',nfrozen

d=(maxdr+1.e-5-mindr)/10.
for j in range(len(drr)):
  a=(drr[j]-mindr)/d
  b=int(a)
  bin[b]=bin[b]+1

print '-------------------------------------'
print 'Distribution displacements:'
for i in range(0,10):
  print mindr+(i+0.5)*d, bin[i]
