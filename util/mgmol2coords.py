# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Python program to generate input coordinates for mgmol from mgmol output
#
# use: python mgmol2coords.py mgmol_output > coords.in
#-------------------------------------------------------------------------------
import sys, string
import outputTools

output=open(sys.argv[1],'r')
Lo    =output.readlines()
lo=len(Lo)  ## no of lines in file

na=outputTools.countNumAtoms(output)
print("#{} atoms".format(na))

searchterm1='IONIC POSITIONS AND FORCES'
searchterm2='Stepper Forces'

previous_coords=[]
coords=[]
names=[]
movables=[]
species=[]
for i in range(0,na):
  previous_coords.append(0)
  coords.append(0)
  movables.append(1)
  species.append(0)

#search for atom names
count=0
list_species=[]
for line in Lo:
  num_matches = line.count('##')
  if num_matches:
    words=line.split()
    if words[0]=='##':
      count=count+1
      if words[1][0]=='*':
        movables[count]=0
        name=words[1][1:]
      else:
        name=words[1]
      if name[1].isdigit():
        sp=name[0]
      else:
        sp=name[0:2]
      names.append(name)
      if sp not in list_species:
        list_species.append(sp)
  if count==na:
    break

   
for line in range(lo): ## loop over lines of file 
  num_matches1 = Lo[line].count(searchterm1)
  num_matches2 = Lo[line].count(searchterm2)
  if num_matches2 or num_matches1:
      j=0
      flag1=0
      flag2=0
      #print 'loop starting at', line+1
      for line2 in range(line+1,lo):
        i=0
        for c in Lo[line2]:
          flag1=0
          if c=='#': 
            if Lo[line2][i+1]=='#':
              flag1=1
              flag2=1 #turn on when first '##' found
              word=Lo[line2][i+3:].split()
              name=word[-7]
              if word[-7][0]=='*':
                name=word[-7][1:]
              x=word[-6]
              y=word[-5]
              z=word[-4]
              #print(name+'\t'+x+'\t'+y+'\t'+z)
              previous_coords[j]=coords[j]
              coords[j]=x.ljust(10)+y.ljust(10)+z.ljust(10)
              j=j+1
              break
          i=i+1
        if flag1!=flag2: break

#for i in range(na):
#  print(coords[i])

#read velocities
searchterm1='IONIC POSITIONS AND VELOCITIES'

found_velocities=0
for line in Lo: ## loop over lines of file 
  num_matches1 = line.count(searchterm1)
  if num_matches1:
    found_velocities=1
    break

vels=[]
if found_velocities:
  for i in range(0,na):
    vels.append("0")

  for line in range(lo): ## loop over lines of file 
    num_matches1 = Lo[line].count(searchterm1)
    if num_matches1:
      j=0
      flag1=0
      flag2=0
      #print 'loop starting at', line+1
      for line2 in range(line+1,lo):
        for c in Lo[line2]:
          flag1=0
          if c=='#': 
            if Lo[line2][1]=='#':
              flag1=1
              flag2=1
              word=Lo[line2].split()
              if len(word)>5:
                vx=word[5]
                vy=word[6]
                vz=word[7]
                #print(name+'\t'+vx+'\t'+vy+'\t'+vz)
                vels[j]=vx.ljust(12)+vy.ljust(12)+vz.ljust(12)
                j=j+1
              break
        if flag1!=flag2: break
else:
  searchterm1='Timestep'
  for line in Lo:
    num_matches1 = line.count(searchterm1)
    if num_matches1:
      words=line.split()
      dt=eval(words[5])
      print("#dt={}".format(dt))
      break  
  
  if previous_coords[0]!=0:
    for i in range(len(coords)):
      words1=coords[i].split()
      words2=previous_coords[i].split()
      vx="%8.6f" % ((eval(words1[0])-eval(words2[0]))/dt)
      vy="%8.6f" % ((eval(words1[1])-eval(words2[1]))/dt)
      vz="%8.6f" % ((eval(words1[2])-eval(words2[2]))/dt)
      vels.append(vx.ljust(11)+vy.ljust(11)+vz.ljust(11))
  

#print new coords file
j=0
for name in names:
  #print(name)
  strname=str(name)
  if len(strname)>1:
    if strname[1].isdigit():
      sp=strname[0]
    else:
      sp=strname[0:2]
  else:
    sp=strname[0]
  #print list_species.index(sp)
  ssp=list_species.index(sp)+1
  print(name+'\t'+str(ssp)+'\t'+coords[j]+'\t'+str(movables[j]),end='')
  if len(vels)>0:
    print('\t'+vels[j])
  else:
    print(' ')
  j=j+1

