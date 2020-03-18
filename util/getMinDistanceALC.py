# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Python program to compare LR centers coordinates in mgmol output
#
# use: python getMinDistanceALC.py mgmol_output
#-------------------------------------------------------------------------------
import sys, string
from math import sqrt
from re import split

input1=open(sys.argv[1],'r')

L1=input1.readlines()
l1=len(L1)  ## no of lines in file


##############################################
# count number of LRs
def getNumLRs(L):
  searchstring='&&'
  found_current_line=0
  already_found_one=0
  na=0
  for line in L: ## loop over lines of file 
    num_matches = string.count(line, searchstring)
    if num_matches:
      found_current_line=1
      already_found_one =1
      na=na+1
    else:
      found_current_line=0
    if found_current_line!=already_found_one:
      break
  return na
##############################################


na1=getNumLRs(L1)

print 'N LRs in file1=', na1


##############################################
def getCoords(coords,L,na):
  l=len(L)  ## no of lines in file
  line_min=0
  searchstring='&&'
  for line in range(l): ## loop over lines of file
    if line>line_min:
      num_matches = string.count(L[line], searchstring)
      if num_matches:
        j=0
        print 'Read coords starting at line ', line
        for line2 in range(line,line+na):
          words=string.split(L[line2])
          x=words[2]
          y=words[3]
          z=words[4]
          coords[j]=x+'\t'+y+'\t'+z
          j=j+1
        line_min=line+na+1
##############################################

# get cell dimension 
cell=[]
searchterm='cell'
for line in L1: ## loop over lines of file 
  num_matches = string.count(line, searchterm)
  if num_matches:
    words=string.split(line)
    words=split("\(",words[0])
    words=split("\)",words[1])
    words=split(",",words[0])
    #print words
    lx=eval(words[0])
    ly=eval(words[1])
    lz=eval(words[2])
 
    print 'cell dimensions=',lx,ly,lz
    break

fac=[0.,-1.,1.]


coords1=[]
for i in range(0,na1):
  coords1.append(0)

getCoords(coords1,L1,na1)

maxdr=0.
minmindr=1000.
drr=[]
bin=[]
for i in range(0,10):
  bin.append(0)

list_distances=[]
na=0
jmin=-1
nfrozen=0
for i in range(na1-1): 
  word1=string.split(coords1[i])
  mindr=1000.
  for j in range(i+1,na1): 
    word2=string.split(coords1[j])
    #check minimage
    dx=10000.
    dy=10000.
    dz=10000.
    for k in range (3):
      dx=min(dx,abs(eval(word1[0])-eval(word2[0])+fac[k]*lx))
    for k in range (3):
      dy=min(dy,abs(eval(word1[1])-eval(word2[1])+fac[k]*ly))
    for k in range (3):
      dz=min(dz,abs(eval(word1[2])-eval(word2[2])+fac[k]*lz))
    dr=sqrt(dx*dx+dy*dy+dz*dz)
    if( dr<mindr ):
      mindr=dr
      jmin=j
    if( mindr<0.001 ):
      break
  
  drr.append(mindr)
  #print mindr, '\t', coords1[i]
      
  list_distances.append(`mindr`+' '+coords1[i]+' '+coords1[jmin])
  if mindr<minmindr:
    minmindr=mindr
    print 'Min. distance so far (for center ',i,' and ',jmin,'):',minmindr
  if mindr>maxdr:
    maxdr=mindr

ang2bohr=1.8897269
bohr2ang=1./ ang2bohr

list_distances.sort()
print '-------------------------------------'
print 'List centers with smallest separating distance:'
#list_distances[na1-10:na1]
for j in range(min(10,na1)):
  word=string.split(list_distances[j])
  print eval(word[0])*bohr2ang,'[Ang]',', center1=',word[1],word[2],word[3],', center2=',word[4],word[5],word[6]

d=(maxdr+1.e-5)/10.
for j in range(len(drr)):
  a=(drr[j])/d
  b=int(a)
  bin[b]=bin[b]+1

