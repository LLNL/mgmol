# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Make list of atoms in xyz format from pdb input 
import sys, string
myfile=open(sys.argv[1],'r')
lines=myfile.readlines()

# list possible species in a dictionary
list_species={'H':0,'He':0,'Li':0,'Be':0,'B':0,'C':0,'N':0,'O':0,'F':0,'Na':0,
              'Mg':0,'Al':0,'Si':0,'P':0,'S':0,'Cl':0,'Ca':0,'Zn':0,'noname':0}
numbers_list=['0','1','2','3','4','5','6','7','8','9']

#count atoms
na=0
for line in lines: ## loop over lines of file 
  if line[0:4] == 'ATOM' or line[0:4] == 'atom' or line[0:6] == 'HETATM':
    na=na+1

minx=1000.
miny=1000.
minz=1000.
maxx=-1000.
maxy=-1000.
maxz=-1000.
for line in lines: ## loop over lines of file 
  if line[0:4] == 'ATOM' or line[0:4] == 'atom' or line[0:6] == 'HETATM':
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    minx=min(x,minx)
    miny=min(y,miny)
    minz=min(z,minz)
    maxx=max(x,maxx)
    maxy=max(y,maxy)
    maxz=max(z,maxz)

print (na)
print ('#generated from file ',sys.argv[1],' ll=(',minx,miny,minz,'),ur=(',maxx,maxy,',',maxz,')')

for line in lines: ## loop over lines of file 
  if line[0:4] == 'ATOM' or line[0:4] == 'atom' or line[0:6] == 'HETATM':
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    name = line[12:16]
    name = name.strip()
    
    words=str.split(line)
      
    firstletter=name[0:1]
    if firstletter in numbers_list:
      firstletter=name[1:2]

    species = firstletter

    if firstletter=='C':
      if name[0:3]=='CLA':
        species="Cl"
      if name[0:2]=='CL':
        species="Cl"
    if firstletter=='S':
      if name[0:3]=='SOD':
        species="Na"
      else:
        if name[0:2]=='Si':
          species="Si"
    if firstletter=='Z':
      if name[0:2]=='ZN':
        species="Zn"

    if species not in list_species.keys():
      print ("ERROR: unknown species", species)
      break

    print (species,'\t',x,'\t',y,'\t',z)
