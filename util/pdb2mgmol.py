# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Make list of atoms in MGmol format from pdb input 
# set movable flag according to occupancy number
import sys, string
myfile=open(sys.argv[1],'r')
lines=myfile.readlines()

# list possible species in a dictionary
list_species={'H':0,'He':0,'Li':0,'Be':0,'B':0,'C':0,'N':0,'O':0,'F':0,'Na':0,
              'Mg':0,'Al':0,'Si':0,'P':0,'S':0,'Cl':0,'Ca':0,'noname':0}
numbers_list=['0','1','2','3','4','5','6','7','8','9']

nel=0
na=0
used_numbers=[]

max_x=-1000.
max_y=-1000.
max_z=-1000.
min_x=1000.
min_y=1000.
min_z=1000.

ang2bohr=1.8897269

nsp=0 #number of species found

for line in lines: ## loop over lines of file 
  if line[0:4] == 'ATOM' or line[0:4] == 'atom' or line[0:6] == 'HETATM':
    x = float(line[30:38])*ang2bohr
    y = float(line[38:46])*ang2bohr
    z = float(line[46:54])*ang2bohr
    name = line[12:16]
    name = name.strip()
    
    occupancy = line[54:60]
    occupancy = occupancy.strip()

    words=string.split(line)
    number  = words[1]
    if number in used_numbers:
      print "ERROR: number ",number," used more than once in pdb input"
      break
    else:
      used_numbers.append(number)
      
    firstletter=name[0:1]
    if firstletter in numbers_list:
      firstletter=name[1:2]

    species = firstletter

    ne=-1
    if firstletter=='H':
      ne=1
    if firstletter=='C':
      if name[0:3]=='CLA':
        species="Cl"
        ne=7
      else:
        ne=4
    if firstletter=='N':
      ne=5
    if firstletter=='O':
      ne=6
    if firstletter=='F':
      ne=7
    if firstletter=='P':
      ne=5
    if firstletter=='S':
      if name[0:3]=='SOD':
        species="Na"
        ne=7
      else:
        if name[0:2]=='Si':
          species="Si"
          ne=4
        else:
          ne=6

    if species in list_species.keys():
      sp=list_species[species]
      if sp==0:
        nsp=nsp+1
        list_species[species]=nsp
        sp=nsp
    else:
      print "ERROR: unknown species", species
      break

    if ne<0:
      print "ERROR: unknown number of electrons for species ", species
      break

    name=species+number

    if x>max_x: max_x=x
    if y>max_y: max_y=y
    if z>max_z: max_z=z
    if x<min_x: min_x=x
    if y<min_y: min_y=y
    if z<min_z: min_z=z

    if ne<0:
      print "could not determine species for atom ", name
      break

    nel=nel+ne
    na=na+1

    print name,'\t',sp,'\t',x,'\t',y,'\t',z,'\t', occupancy[0:1]

print "#Number of atoms     =",na
print "#Number of electrons =",nel
print "#Max. x=", max_x
print "#Max. y=", max_y
print "#Max. z=", max_z
print "#Min. x=", min_x
print "#Min. y=", min_y
print "#Min. z=", min_z
