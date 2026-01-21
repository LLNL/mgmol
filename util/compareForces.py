# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Python program to compare forces between 2 mgmol outputs
#
# use: python diff_forces.py mgmol_output1 mgmol_output2
#-------------------------------------------------------------------------------
import sys, string
from math import sqrt
import matplotlib.pyplot as plt

input1=open(sys.argv[1],'r')
input2=open(sys.argv[2],'r')

frame=-1
if len(sys.argv)>3:
  frame=eval(sys.argv[3])
  print( 'Input argument: Frame=',frame )

lines1=input1.readlines()
lines2=input2.readlines()

##############################################
# count number atoms
def getNumAtoms(lines):
  found_current_line=0
  already_found_one=0
  na=0
  flag=0
  for line in lines: ## loop over lines of file
    if line.count('FORCES') or line.count('Forces'):
      flag=1
    if line.count('## ') & flag==1:
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


na1=getNumAtoms(lines1)
na2=getNumAtoms(lines2)

print( 'N atoms in file1=', na1)
print( 'N atoms in file2=', na2)

##############################################

def getForces(names,coords,forces,lines,fframe):
  na=len(names)
  l=len(lines)  ## no of lines in file
  line_min=0  
  cur_frame=0
  for line in range(l): ## loop over lines of file1 
    if line>line_min:
      thisline = lines[line]
      num_matches3 = thisline.count('md_ProjectForces')
      if num_matches3:
        return

      num_matches1 = thisline.count('Forces')
      num_matches2 = thisline.count('FORCES')
      num_matches3 = thisline.count('Timer:')
      num_matches=num_matches1+num_matches2
      if num_matches and not num_matches3:
        j=0
        print( 'Frame: ',cur_frame,' Read forces starting at line ', line,' for ',na,' atoms')
        for line2 in range(line+1,line+na+2):
          thisline2 = lines[line2]
          if thisline2.count('##')>0:
            words=thisline2.split()
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
          print ('break')
          break
        cur_frame=cur_frame+1

##############################################

forces1=[0]*na1
coords1=[0]*na1
names1=[0]*na1

forces2=[0]*na2
coords2=[0]*na2
names2=[0]*na2
  
getForces(names1,coords1,forces1,lines1,frame)
getForces(names2,coords2,forces2,lines2,frame)

mindf=100.
maxdf=0.
avg=0.
avgx=0.
avgy=0.
avgz=0.
imax=0
jmax=0
dff=[]
bins=[0] * 10

na=0
for i in range(na1):
  force1=forces1[i]
  word1=force1.split()
  for j in range(na2):
    force2=forces2[j]
    word2=force2.split()
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
      print (names1[i],': delta f=',df)

print ('na=',na)
avg=avg/na
avgx=avgx/na
avgy=avgz/na
avgz=avgz/na

print ('N atoms = ', na)
print ('Avg. df = ',avgx,avgy,avgz)
print ('Avg. |df| = ',avg)
print ('Min. |df| = ',mindf)
print ('Max. |df| = ',maxdf)
print ('Atoms with largest force difference:')
filename1=sys.argv[1]
filename1=filename1.ljust(15)
filename2=sys.argv[2]
filename2=filename2.ljust(15)
print (filename1,'\t',names1[imax],'\t',coords1[imax],'\t',forces1[imax])
print (filename2,'\t',names2[jmax],'\t',coords2[jmax],'\t',forces2[jmax])

delf=(maxdf+1.e-5-mindf)/10.
for j in range(na):
  a=(dff[j]-mindf)/delf
  b=int(a)
  bins[b]=bins[b]+1

for i in range(0,10):
  print (mindf+(i+0.5)*delf, bins[i])

plt.hist(dff, bins=10, edgecolor="black")
plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.xlabel('force error magnitude [Ha/Bohr]',fontsize=12)
plt.ylabel('frequency',fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('errorForces.png', dpi=100)
