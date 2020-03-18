# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
import sys, string

# list possible species in a dictionary with their count
species={'H':0,'D':0,'He':0,'Li':0,'Be':0,'B':0,'C':0,'N':0,'O':0,'F':0,'Na':0,
         'Mg':0,'Al':0,'Si':0,'P':0,'S':0,'Cl':0,'Ca':0,'In':0,'Au':0}

ang2bohr=1.8897269
bohr2ang=1./ang2bohr

########################################################################
def makePDBobject(na,anames,acoords,aocc,pdbobject):
  for i in range(na):
    pdbobject.append('ATOM')

  #atom number in columns 7-11
  for i in range(na):
    #number=str(i)
    if anames[i][0:2] in species.keys():
      number=anames[i][2:]
    else:
      number=anames[i][1:]
    number=number.rjust(5)
    pdbobject[i]=pdbobject[i]+'  '+number

  #atom name in columns 13-16
  for i in range(na):
    if anames[i][0:2] in species.keys():
      name=anames[i][0:2]
    else:
      name=anames[i][0:1]
    name=name.ljust(4)
    pdbobject[i]=pdbobject[i]+' '+name

  #atom name in columns 18-20
  resname='   '
  for i in range(na):
    pdbobject[i]=pdbobject[i]+' '+resname

  #atom number in columns 23-26
  resnum='    '
  for i in range(na):
    pdbobject[i]=pdbobject[i]+'  '+resnum

  #x in columns 31-38  
  #y in columns 39-46 
  #z in columns 47-54
  for i in range(na): 
    word=string.split(acoords[i])
    pdbobject[i]=pdbobject[i]+'    '
    for j in range(3):
      x=eval(word[j])
      x = "%.3f" % (x*bohr2ang)
      x=x.rjust(8)
      pdbobject[i]=pdbobject[i]+x
    occ=str(aocc[i]).rjust(4)
    pdbobject[i]=pdbobject[i]+occ
  
########################################################################
def printPDBfile(natoms,anames,acoords,aocc,filename):
  lineout=[]
  makePDBobject(natoms,anames,acoords,aocc,lineout)

  if filename:
    output = open(filename,'w')
    for i in range(0,natoms):
      output.write(lineout[i])
      output.write('\n')
  else:
    for line in lineout:
      print line
