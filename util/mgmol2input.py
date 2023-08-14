# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Python program to generate input files for mgmol from mgmol output
#
# use: python mgmol2input.py mgmol_input mgmol_output > mgmol_new_input
#-------------------------------------------------------------------------------
import sys, string
import outputTools

ifile =open(sys.argv[1],'r')
output=open(sys.argv[2],'r')
Lo    =output.readlines()
lo=len(Lo)  ## no of lines in file

# list possible species in a dictionary with their count
species={'H':0,'He':0,'Li':0,'Be':0,'B':0,'C':0,'N':0,'O':0,'F':0,'Na':0,
         'Mg':0,'Al':0,'Si':0,'P':0,'S':0,'Cl':0,'Ca':0,'In':0,'Au':0}

##############################################
# count number orbitals
def getNumOrbitals(lines):
  for line in lines: ## loop over lines of file
    num_matches = line.count('Number of states')
    if num_matches:
      word=line.split()
      no=eval(word[4])
      break
  return no
##############################################

na=outputTools.countNumAtoms(output)
print("#{} atoms".format(na))

searchterm1='IONIC POSITIONS AND FORCES'
searchterm2='Stepper Forces'

previous_coords=[]
coords=[]
names=[]
for i in range(0,na):
  previous_coords.append(0)
  coords.append(0)
  names.append(0)
  
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
              flag2=1
              word=Lo[line2][i+3:].split()
              name=word[-7]
              if word[-7][0]=='*':
                name=word[-7][1:]
              x=word[-6]
              y=word[-5]
              z=word[-4]
              #print name+'\t'+x+'\t'+y+'\t'+z
              names[j]=name
              previous_coords[j]=coords[j]
              coords[j]=x.ljust(10)+y.ljust(10)+z.ljust(10)
              j=j+1
              break
          i=i+1
        if flag1!=flag2: break

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
                #print name+'\t'+vx+'\t'+vy+'\t'+vz
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
  

#read MLWF centers
searchterm1='centers'
searchterm2='spreads'
print("#{} names".format(len(names)))

no=getNumOrbitals(Lo)
print("#{} orbitals".format(no))
lrs=[]
for i in range(0,no):
  lrs.append(0)
  
for line in range(lo): ## loop over lines of file 
  num_matches1 = Lo[line].count(searchterm1)
  if num_matches1:
    num_matches2 = Lo[line].count(searchterm2)
    if num_matches2:
      j=0
      found_current_line=0
      already_found_one =0
      #print 'loop starting at', line+1
      for line2 in range(line+1,line+no+5):
        words=Lo[line2].split()
        if len(words)>0:
          if words[0]=='&&':
            found_current_line=1
            already_found_one =1
            x=words[2]
            y=words[3]
            z=words[4]
            s=words[5]
            lrs[j]=x.ljust(12)+y.ljust(12)+z.ljust(12)
            j=j+1
          else:
            found_current_line=0
        if found_current_line!=already_found_one:
          break

#print new input file
lines=ifile.readlines()
natoms_printed=0
count_orbital=0
for line in lines: ## loop over lines of file
  word=line.split()
  new_line=0
  if natoms_printed<na:
    if len(word)>0:
      name=word[0]
      #print(name)
      for i in range(0,na):
        if name==names[i]:
          movable=1
          if len(word)>5:
            movable=word[5]
          print(name+'\t'+word[1]+'\t'+coords[i]+'\t'+str(movable),end='')
          if len(vels)>0:
            print('\t'+vels[i])
          else:
            print(" ")
          natoms_printed=natoms_printed+1
          new_line=1
          break
  
  else:
    if count_orbital<no:
      if len(word)>0:
        if word[0][0]!='#' and word[0][0]!=' ':
          if len(word)>=3:
            if len(word)>3:
              print(lrs[count_orbital],word[3])
            else:
              print(lrs[count_orbital])
            count_orbital=count_orbital+1
            new_line=1
  
  if new_line==0:
    print(line)

if count_orbital<no:
  print("#LRs {}".format(count_orbital))
  for jj in range(count_orbital,no):
    print(lrs[jj])

