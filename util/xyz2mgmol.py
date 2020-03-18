# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Make list of atoms in MGmol format from xyz input 
import sys, string
myfile=open(sys.argv[1],'r')
lines=myfile.readlines()

# list possible species in a dictionary
list_species={'H':0,'He':0,'Li':0,'Be':0,'B':0,'C':0,'N':0,'O':0,'F':0,'Na':0,
              'Mg':0,'Al':0,'Si':0,'P':0,'S':0,'Cl':0,'Ca':0,'Zn':0,'Mo':0,
              'noname':0}

charges={'H':1,'Li':3,'Be':4,'B':3,'C':4,'N':5,'O':6,'F':7,'Na':7,
         'Mg':-1,'Al':3,'Si':4,'P':5,'S':6,'Cl':7,'Zn':12,'Mo':14}

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

count=0
number=0
for line in lines: ## loop over lines of file 
  count=count+1
  if count>2:
    number=number+1
    words=line.split()

    species = words[0]
    if words[0]=='SOD':
        species="Na"

    ne=charges[species]
    if ne<0:
      print("ERROR: unknown number of electrons for species {}".format(species))
      break

    if species in list_species.keys():
      sp=list_species[species]
      if sp==0:
        nsp=nsp+1
        list_species[species]=nsp
        sp=nsp
    else:
      print("ERROR: unknown species {}".format(species))
      break

    if ne<0:
      print("ERROR: unknown number of electrons for species {}".format(species))
      break

    x = eval(words[1])*ang2bohr
    y = eval(words[2])*ang2bohr
    z = eval(words[3])*ang2bohr

    name=species+str(number)

    if x>max_x: max_x=x
    if y>max_y: max_y=y
    if z>max_z: max_z=z
    if x<min_x: min_x=x
    if y<min_y: min_y=y
    if z<min_z: min_z=z

    if ne<0:
      print("could not determine species for atom {}".format(name))
      break

    nel=nel+ne
    na=na+1

    print(name.ljust(7), str(sp).rjust(3), \
          str(x).rjust(16),str(y).rjust(16),str(z).rjust(16))

print("#Number of atoms     = {}".format(na))
print("#Number of electrons = {}".format(nel))
print("#Max. x={}".format(max_x))
print("#Max. y={}".format(max_y))
print("#Max. z={}".format(max_z))
print("#Min. x={}".format(min_x))
print("#Min. y={}".format(min_y))
print("#Min. z={}".format(min_z))
print("#lx[Bohr]={}".format(max_x-min_x))
print("#ly[Bohr]={}".format(max_y-min_y))
print("#lz[Bohr]={}".format(max_z-min_z))

