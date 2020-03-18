# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Make list of LR centers in MGmol format from xyz input 
import sys, string
from math import sqrt

myfile=open(sys.argv[1],'r')
lines=myfile.readlines()

# list possible species in a dictionary
#ele_nb={'H':1,'Li':1,'C':4,'N':3,'O':6,'F':7,'Si':4,'P':5,'S':6,'noname':-1}
#extra e- on P, compensated by Li, assuming same number of Li and P
ele_nb={'H':1,'Li':0,'C':4,'N':3,'O':6,'F':7,'Si':4,'P':6,'S':6,'noname':-1}

# covalent radius (from RasMol, abstree.h, first number/250)
cov_radii={'H':0.32,'Li':0.68,'C':0.72,'N':0.68,'O':0.68,'F':0.64,'Na':0.972,
           'Mg':1.1,'Si':1.2,
           'P':1.036,'S':1.02,'Cl':1.,'K':1.328,'Ca':0.992,'noname':-1.}
numbers_list=['0','1','2','3','4','5','6','7','8','9']


lx=0.
ly=0.
lz=0.
maxd=0.

ang2bohr=1.8897269

atoms_ne={}
atoms_coords={}
pairs=[]

########################################
def get_sp(name):
  sp=name[0:2]
  if sp[1] in numbers_list:
    sp=name[0]

  return sp

########################################
def distance(x1,y1,z1,x2,y2,z2):
  mind=10000.
  dx=(x1-x2)
  dy=(y1-y2)
  dz=(z1-z2)
  d=sqrt(dx*dx+dy*dy+dz*dz)
  if d<maxd:
    return d
  
  ab=[0,1,-1]
  for i in ab:
    for j in ab:
      for k in ab:
        dx=(x1-x2)+i*lx
        dy=(y1-y2)+j*ly
        dz=(z1-z2)+k*lz
        d=sqrt(dx*dx+dy*dy+dz*dz)
        mind=min(mind,d)
        if mind<maxd:
          return mind
  return mind

########################################
def adjacent(name1,name2):
  d=distancePair(name1,name2)
  
  sp1=get_sp(name1)
  sp2=get_sp(name2)

  #print sp1, sp2, d
  if d<(cov_radii[sp1]+cov_radii[sp2]+0.56):
    return True
    
  return False

########################################
def distancePair(name1,name2):
  words1=string.split(atoms_coords[name1])
  words2=string.split(atoms_coords[name2])
  x1=eval(words1[0])
  y1=eval(words1[1])
  z1=eval(words1[2])
  x2=eval(words2[0])
  y2=eval(words2[1])
  z2=eval(words2[2])
  
  return distance(x1,y1,z1,x2,y2,z2)

########################################
def assign_pair(name,rname,maxcount):
  if name<rname:
    pair=name+' '+rname
  else:
    pair=rname+' '+name
  
  if pairs.count(pair)<maxcount:
    print rname, name,
    words1=string.split(atoms_coords[name])
    words2=string.split(atoms_coords[rname])
    x1=eval(words1[0])
    y1=eval(words1[1])
    z1=eval(words1[2])
    x2=eval(words2[0])
    y2=eval(words2[1])
    z2=eval(words2[2])
    d=distance(x1,y1,z1,x2,y2,z2)
    x=0.5*(x1+x2)
    y=0.5*(y1+y2)
    z=0.5*(z1+z2)
    print x,y,z,d
    pairs.append(pair)
    atoms_ne[name]=atoms_ne[name]-1
    atoms_ne[rname]=atoms_ne[rname]-1

########################################
def doubleBonds(atoms_list,bonded_atoms):
  search_list=list(set(atoms_list))
  i=0
  nb=len(bonded_atoms)
  na=-1
  while na<len(atoms_list):
    na=len(atoms_list)
    print '#nested=',i
    #print len(atoms_list),len(bonded_atoms)
    
    last_list=[]
    #print 'search list: ', search_list
    for name in search_list:
      #print name
      
      #build list of adjacent atoms
      adjacent_list=[]
      for rname in bonded_atoms:
        if rname!=name:
          if adjacent(name,rname):
            adjacent_list.append(rname)
            last_list.append(rname)
            atoms_list.append(rname)
      
      if len(adjacent_list)>0:
        #print 'name=',name
        #print adjacent_list
        if atoms_ne[name]>0:
          dmin=10000.
          closest_name=""
          for rname in adjacent_list:
            if atoms_ne[rname]>0:
              d=distancePair(name,rname)
              print 'name,rname,d',name,rname,d
              if d<dmin:
                dmin=d
                closest_name=rname
          #print 'closest:',closest_name, dmin
          if dmin<10.:
            assign_pair(name,closest_name,2)
    
    search_list=list(set(last_list))
    atoms_list=list(set(atoms_list))
    
    print '#Num. atoms: ',len(atoms_list)
    i=i+1

  bonded_atoms=list(set(bonded_atoms))
  return
############################

number=0

words=string.split(lines[1])
lx=eval(words[0])
ly=eval(words[1])
lz=eval(words[2])
print '#Domain: ',lx,ly,lz
maxd=min(0.5*lx,0.5*ly)
maxd=min(maxd,0.5*lz)
#maxd=0.5*sqrt(lx*lx+ly*ly+lz*lz)

## loop over lines of file
for line in lines[2:]:
  #print line,
  words=string.split(line)
  
  if len(words)<1:
    break

  number=number+1
  species = words[0]

  ra=cov_radii[species]
  ne=ele_nb[species]
  if ne<0:
    print "ERROR: unknown number of electrons for species ", species
    break
  if ra<0:
    print "ERROR: unknown covariant radius for species ", species
    break

  x = words[1]
  y = words[2]
  z = words[3]
  
  name=species+`number`
  atoms_coords[name]=(x+'\t'+y+'\t'+z)
  atoms_ne[name]=ne
  #print name

bonded_atoms=[]

#first look for atoms bonds with H involved
bonded_atomsH=[]
for name in atoms_ne.keys():  
  sp=get_sp(name)
  if sp=='H':
    for rname in atoms_ne.keys():
      if rname!=name:
        if adjacent(name,rname):
          if atoms_ne[name]>0 and atoms_ne[rname]>0:
            assign_pair(name,rname,1)
            bonded_atoms.append(name)
            bonded_atoms.append(rname)
            bonded_atomsH.append(rname)

bonded_atomsH=list(set(bonded_atomsH))
            
#first look for atoms bonds with O involved
bonded_atomsO=[]
for name in atoms_ne.keys():  
  sp=get_sp(name)
  if sp=='O':
    for rname in atoms_ne.keys():
      if rname!=name:
        if adjacent(name,rname):
          if atoms_ne[name]>0 and atoms_ne[rname]>0:
            assign_pair(name,rname,1)
            bonded_atoms.append(name)
            bonded_atoms.append(rname)
            bonded_atomsO.append(rname)

bonded_atomsO=list(set(bonded_atomsO))
            
bonded_atoms=list(set(bonded_atoms))
#print 'List H:',bonded_atomsH

#look for single bonds
for name in atoms_ne.keys():
  if atoms_ne[name]>0:
    for rname in atoms_ne.keys():
      if rname!=name:
        if adjacent(name,rname):
          if atoms_ne[rname]>0:
            assign_pair(name,rname,1)
            bonded_atoms.append(name)
            bonded_atoms.append(rname)

bonded_atoms=list(set(bonded_atoms))
#print 'List single bonded atoms:',bonded_atoms
print '#Number of single bonded atoms:',len(bonded_atoms)

#put pairs of electrons on atoms
print '#Pairs on atoms '
for i in range(3):
  for name in atoms_ne.keys():
    if atoms_ne[name]>1:
      #print 'Pair on atom ',name
      print name,atoms_coords[name]
      atoms_ne[name]=atoms_ne[name]-2

#look for second bond for "O"
print '#Double bonds for O'
for name in atoms_ne.keys():  
  sp=get_sp(name)
  if sp=='O' and atoms_ne[name]>0:
    for rname in atoms_ne.keys():
      if rname!=name and atoms_ne[rname]>0:
        if adjacent(name,rname):
          assign_pair(name,rname,2)
          bonded_atoms.append(name)
          bonded_atoms.append(rname)

bonded_atoms=list(set(bonded_atoms))

#look for double bonds
print '#Look for double bonds starting from O..'
atoms_list=bonded_atomsO
doubleBonds(atoms_list,bonded_atoms)

print '#Look for double bonds again, starting from H..'
atoms_list=bonded_atomsH
doubleBonds(atoms_list,bonded_atoms)


#look for remaining bonds
for name in atoms_ne.keys():
  if atoms_ne[name]>0:
    for rname in atoms_ne.keys():
      if rname!=name:
        if adjacent(name,rname):
          if atoms_ne[rname]>0:
            assign_pair(name,rname,1)
            bonded_atoms.append(name)
            bonded_atoms.append(rname)

print '#Unpaired electrons:'
for name in atoms_ne.keys():
  if atoms_ne[name]>0:
    print name,atoms_ne[name]

#print atoms_ne

pairs_dic={}
print '###Pairs###'
for pair in pairs:
  pairs_dic[pair]=pairs.count(pair)

for pair in pairs_dic.keys():
  print pair, pairs_dic[pair]

print "#Number of atoms     =",number

