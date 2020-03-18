# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# use: python getDistributionOrbSpreads.py mgmol_output
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

def getSpreads(spreads,L,na):
  searchstring_start='ABPG'
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
        flag=0
        j=0
        print 'Read ',na,' spreads starting at line ', l
        for line2 in range(l,l+na):
          words=string.split(L[line2])
          s=words[5]
          spreads[j]=s
          j=j+1
        line_min=l+na+1
    l=l+1
##############################################


na1=getNumLRs(L1)

print 'N LRs in file=', na1


##############################################

spreads1=[]
for i in range(0,na1):
  spreads1.append(0)

getSpreads(spreads1,L1,na1)

list_spreads=[]
for i in range(na1): 
  word1=string.split(spreads1[i])  
  list_spreads.append(word1[0])

##############################################
ang2bohr=1.8897269
bohr2ang=1./ ang2bohr

list_spreads.sort()
print '-------------------------------------'
mins=eval(list_spreads[0])
maxs=eval(list_spreads[na1-1])
print 'Min spread = ',mins
print 'Max spread = ',maxs

nbins=10
dx=(maxs+0.01-mins)/nbins

dist=[]
for i in range(nbins):
  dist.append(0)
for j in range(na1):
  n=int((eval(string.split(list_spreads[j])[0])-mins)/dx)
  dist[n]=dist[n]+1

for i in range(nbins):
  print mins+dx*(i+0.5),'\t',dist[i]
