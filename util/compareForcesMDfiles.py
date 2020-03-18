# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Python program to compare forces between 2 mgmol outputs
#
# use: python compareForcesMDfiles.py output1 output2
#-------------------------------------------------------------------------------
import sys, string
from math import sqrt

input1=open(sys.argv[1],'r')
input2=open(sys.argv[2],'r')

lines1=input1.readlines()
lines2=input2.readlines()

##############################################
# count number atoms
def getNumAtoms(lines):
  na=0
  for line in lines: ## loop over lines of file 
    words=string.split(line)
    if len(words)==10:
      na=na+1

  return na
##############################################


na1=getNumAtoms(lines1)
na2=getNumAtoms(lines2)

print 'N atoms in file1=', na1
print 'N atoms in file2=', na2


##############################################

def getForces(names,coords,forces,lines):
  j=0
  for line in lines: ## loop over lines of file1 
    words=string.split(line)
    if len(words)==10:

      name=words[0]
      #print name
      x=words[1]
      y=words[2]
      z=words[3]
      fx=words[4]
      fy=words[5]
      fz=words[6]
      names[j]=name
      coords[j]=x+'\t'+y+'\t'+z
      forces[j]=fx+'\t'+fy+'\t'+fz
      j=j+1

##############################################

forces1=[]
coords1=[]
names1=[]
for i in range(0,na1):
  forces1.append(0)
  coords1.append(0)
  names1.append(0)

forces2=[]
coords2=[]
names2=[]
for i in range(0,na2):
  forces2.append(0)
  coords2.append(0)
  names2.append(0)
  
  
getForces(names1,coords1,forces1,lines1)
getForces(names2,coords2,forces2,lines2)

mindf=100.
maxdf=0.
avg=0.
avgx=0.
avgy=0.
avgz=0.
imax=0
jmax=0
dff=[]
bin=[]
for i in range(0,10):
  bin.append(0)

##############################################
def subtractAverageForce(forces):
  avgx=0.
  avgy=0.
  avgz=0.
  na=len(forces)
  for i in range(na): 
    word=string.split(forces[i])
    fx=eval(word[0])
    fy=eval(word[1])
    fz=eval(word[2])
    avgx=avgx+fx
    avgy=avgy+fy
    avgz=avgz+fz

  avgx=avgx/na
  avgy=avgy/na
  avgz=avgz/na
  
  for i in range(na):
    word=string.split(forces[i])
    fx=eval(word[0])-avgx
    fy=eval(word[1])-avgy
    fz=eval(word[2])-avgz
    forces[i]=`fx`+'\t'+`fy`+'\t'+`fz`

##############################################

#subtractAverageForce(forces1)
#subtractAverageForce(forces2)

na=0
for i in range(na1): 
  word1=string.split(forces1[i])
  for j in range(na2): 
    word2=string.split(forces2[j])
    if names1[i]==names2[j]:
      fx1=eval(word1[0])
      fy1=eval(word1[1])
      fz1=eval(word1[2])
      f1=sqrt(fx1*fx1+fy1*fy1+fz1*fz1)
      fx2=eval(word2[0])
      fy2=eval(word2[1])
      fz2=eval(word2[2])
      dfx=fx1-fx2
      dfy=fy1-fy2
      dfz=fz1-fz2
      df=sqrt(dfx*dfx+dfy*dfy+dfz*dfz)
      dff.append(df)
      avg=avg+df
      avgx=avgx+dfx
      avgy=avgy+dfy
      avgz=avgz+dfz
      if df>maxdf:
        maxdf=df
        imax=i
        jmax=j
      if df<mindf:
        mindf=df
      na=na+1
      print names1[i],': delta f=',df
print 'na=',na      
avg=avg/na

print 'N atoms =', na
print 'Avg. df=',avgx,avgy,avgz
print 'Avg. |df|=',avg
print 'Min. df=',mindf
print 'Max. df=',maxdf
print 'df max for atom ',names1[imax],' and ',names2[jmax]
print 'Forces atoms with largest force difference:'
filename1=sys.argv[1]
filename1=filename1.ljust(15)
filename2=sys.argv[2]
filename2=filename2.ljust(15)
print filename1,'\t',names1[imax],'\t',coords1[imax],'\t',forces1[imax]
print filename2,'\t',names2[jmax],'\t',coords2[jmax],'\t',forces2[jmax]

delf=(maxdf+1.e-5-mindf)/10.
for j in range(na):
  a=(dff[j]-mindf)/delf
  b=int(a)
  bin[b]=bin[b]+1

for i in range(0,10):
  print mindf+(i+0.5)*delf, bin[i]

#for j in range(na): 
#  print dff[j]
