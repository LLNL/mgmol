# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Python program to measure distances between localization centers
# in MGmol output and detect centers close to each other
#
# use: python minDistanceALC.py mgmol_output
#-------------------------------------------------------------------------------
import sys, string
from math import sqrt

input=open(sys.argv[1],'r')

L1=input.readlines()
l1=len(L1)  ## no of lines in file


##############################################
# count number of LRs
def getNumLRs(L):
  searchstring_start='%%'
  searchstring='&&'
  found_current_line=0
  already_found_one=0
  na=0
  flag=0
  for line in L: ## loop over lines of file 
    num_matches_start = string.count(line, searchstring_start)
    if num_matches_start:
      flag=1
    if flag>0:
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

def getCoords(coords,L,na):
  searchstring_start='%%'
  line_min=0
  searchstring='&&'
  flag=0
  l=0
  for line in L: ## loop over lines of file
    num_matches_start = string.count(line, searchstring_start)
    if num_matches_start:
      flag=1
    if flag>0 and l>line_min:
      num_matches = string.count(line, searchstring)
      if num_matches:
        j=0
        print 'Read ',na,' coords starting at line ', l
        for line2 in range(l,l+na):
          words=string.split(L[line2])
          x=words[2]
          y=words[3]
          z=words[4]
          coords[j]=x+'\t'+y+'\t'+z
          j=j+1
        line_min=l+na+1
    l=l+1
##############################################


na1=getNumLRs(L1)

print 'N LRs in file=', na1


##############################################
# get cell dimension 
cell=[]
searchterm='Dimension'
for line in L1: ## loop over lines of file 
  num_matches = string.count(line, searchterm)
  if num_matches:
    word=string.split(line)
 
    lx=eval(word[2])
    ly=eval(word[3])
    lz=eval(word[4])
    
    break
 
print 'cell dimensions=',lx,ly,lz

##############################################

coords1=[]
for i in range(0,na1):
  coords1.append(0)

getCoords(coords1,L1,na1)

minmindr=1000.
drr=[]
bin=[]
for i in range(0,10):
  bin.append(0)

list_distances=[]
fac=[0.,-1.,1.]
na=0
jmin=-1
nfrozen=0
for i in range(1,na1): 
  word1=string.split(coords1[i])
  mindr=1000.
  for j in range(i): 
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
      word2_save=word2
    if( mindr<0.001 ):
      break
  
  drr.append(mindr)
  #print mindr, '\t', coords1[i]
      
  list_distances.append(`mindr`+' '+coords1[i]+' '+coords1[jmin])
  if mindr<minmindr:
    minmindr=mindr
    print 'Min. distance so far:',minmindr
    print 'center1=',word1,     ', i=',i
    print 'center2=',word2_save,', j=',j

##############################################
ang2bohr=1.8897269
bohr2ang=1./ ang2bohr

list_distances.sort()
print '-------------------------------------'
print 'List centers with smallest separating distance:'
for j in range(min(10,na1)):
  word=string.split(list_distances[j])
  print eval(word[0])*bohr2ang,'[Ang]',
  print ', center1=',word[1],word[2],word[3],
  print ', center2=',word[4],word[5],word[6]
