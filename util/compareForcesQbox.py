# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Python program to compare forces between 1 mgmol output and a Qbox output
#
# use: python diff_forces.py mgmol_output qbox_output
#-------------------------------------------------------------------------------
import sys, string
from math import sqrt

input1=open(sys.argv[1],'r')
inputQbox=open(sys.argv[2],'r')

frame=-1
if len(sys.argv)>3:
  frame=eval(sys.argv[3])
  print 'Input argument: Frame=',frame

L1=input1.readlines()
L2=inputQbox.readlines()

star='*'

##############################################
# count number atoms
def getNumAtoms(L):
  searchterm1='## '
  searchterm2='FORCES'
  searchterm3='Forces'
  found_current_line=0
  already_found_one=0
  na=0
  flag=0
  for line in L: ## loop over lines of file 
    num_matches1 = string.count(line, searchterm1)
    num_matches2 = string.count(line, searchterm2)
    num_matches3 = string.count(line, searchterm3)
    if num_matches2 or num_matches3:
      flag=1
    if num_matches1 & flag==1:
      #print 'line=',line
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

print 'N atoms in file1=', na1


##############################################

searchterm1='Forces'
searchterm2='FORCES'
searchterm3='md_ProjectForces'

def getForces(names,coords,forces,L,fframe):
  na=len(names)
  l=len(L)  ## no of lines in file
  line_min=0  
  cur_frame=0
  for line in range(l): ## loop over lines of file1 
    if line>line_min:
      num_matches3 = string.count(L[line], searchterm3)
      if num_matches3:
        return

      num_matches1 = string.count(L[line], searchterm1)
      num_matches2 = string.count(L[line], searchterm2)
      num_matches3 = string.count(L[line], 'Timer:')
      num_matches=num_matches1+num_matches2
      if num_matches and not num_matches3:
        j=0
        print 'Frame: ',cur_frame,' Read forces starting at line ', line,' for ',na,' atoms'
        for line2 in range(line+1,line+na+2):
          if string.count(L[line2], '##')>0:
            words=string.split(L[line2])
            shift=0
            while words[shift]!='##':
              shift=shift+1
            shift=shift+1
            if words[shift]=='*':
              shift=shift+1
            word=words[shift:]
            name=word[0]
            #print name
            x=word[1]
            y=word[2]
            z=word[3]
            fx=word[4]
            fy=word[5]
            fz=word[6]
            names[j]=name
            coords[j]=x+'\t'+y+'\t'+z
            forces[j]=fx+'\t'+fy+'\t'+fz
            j=j+1
        line_min=line+na+2
        if fframe==cur_frame:
          print 'break'
          break
        cur_frame=cur_frame+1

##############################################

def getForcesQbox(names,coords,forces,L):
  
  j=0
  name=""
  for line in L: ## loop over lines of file 
      num_matches = string.count(line, "<atomset>")
      if num_matches:
        j=0
      num_matches = string.count(line, "<atom name")
      if num_matches:
         words=string.split(line)
         #start = string.find('"')
         name=words[1][6:-1]
         names[j]=name
      num_matches = string.count(line, "<position>")
      if num_matches:
         words=string.split(line)
         x=words[1]
         y=words[2]
         z=words[3]
         coords[j]=x+'\t'+y+'\t'+z
      
      num_matches = string.count(line, "<force>")
      if num_matches:
         words=string.split(line)
         fx=words[1]
         fy=words[2]
         fz=words[3]
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
for i in range(0,na1):
  forces2.append(0)
  coords2.append("")
  names2.append("")
  
  
getForces(names1,coords1,forces1,L1,frame)
getForcesQbox(names2,coords2,forces2,L2)

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

na=0
for i in range(na1): 
  word1=string.split(forces1[i])
  fx1=eval(word1[0])
  fy1=eval(word1[1])
  fz1=eval(word1[2])
  f1=sqrt(fx1*fx1+fy1*fy1+fz1*fz1)
  for j in range(na1): 
    if names2[j]==names1[i]:
      word2=string.split(forces2[j])
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
print 'df max for atom ',names1[imax]
print 'Forces atoms with largest force difference:'
filename1=sys.argv[1]
filename1=filename1.ljust(15)
filename2=sys.argv[2]
filename2=filename2.ljust(15)
print filename1,'\t',names1[imax],'\t',coords1[imax],'\t',forces1[imax]
print filename2,'\t',names2[jmax],'\t',coords2[jmax],'\t',forces2[imax]

delf=(maxdf+1.e-5-mindf)/10.
for j in range(na):
  a=(dff[j]-mindf)/delf
  b=int(a)
  bin[b]=bin[b]+1

for i in range(0,10):
  print mindf+(i+0.5)*delf, bin[i]

#for j in range(na): 
#  print dff[j]
