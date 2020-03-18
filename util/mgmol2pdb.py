# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Python program to generate pdb files from mgmol output
#
# use:
# python mgmol2pdb.py mgmol_output > pdb
# python mgmol2pdb.py mgmol_output dump_freq
#-------------------------------------------------------------------------------
import sys, string
from numpy import array
import outputTools

input_   =open(sys.argv[1],'r')
filename_=sys.argv[1]

default_dump_freq_=100000
dump_freq_=default_dump_freq_ # interval between frames for files generation
if len(sys.argv)>2:
  dump_freq_=eval(sys.argv[2])

L=input_.readlines()
l=len(L)  ## no of lines in file

# list possible species in a dictionary with their count
species={'H':0,'D':0,'He':0,'Li':0,'Be':0,'B':0,'C':0,'N':0,'O':0,'F':0,'Na':0,
         'Mg':0,'Al':0,'Si':0,'P':0,'S':0,'Cl':0,'Ca':0,'In':0,'Au':0}

ang2bohr=1.8897269
bohr2ang=1./ang2bohr

na=outputTools.countNumAtoms(input_)
#print 'N atoms=', na

searchterm='Origin'
ox=0.
oy=0.
oz=0.
for line in L:
  num_matches = string.count(line, searchterm)
  if num_matches:
    words=string.split(line)
    ox=eval(words[2])
    oy=eval(words[3])
    oz=eval(words[4])
origin=array([ox,oy,oz])

searchterm='Dimension'
lx=0.
ly=0.
lz=0.
for line in L:
  num_matches = string.count(line, searchterm)
  if num_matches:
    words=string.split(line)   
    lx=eval(words[2])
    ly=eval(words[3])
    lz=eval(words[4])

cell=array([lx,ly,lz])
end=array([lx,ly,lz])

########################################################################
def readAtomicPositions(first_line,last_line,anames,acoords):
  j=0
  flag1=0
  flag2=0
  #print 'loop starting at', first_line+1
  for line2 in range(first_line,last_line):
    i=0
    for c in L[line2]:
      flag1=0
      if c=='#': 
        if L[line2][i+1]=='#':
          flag1=1
          flag2=1
          word=string.split(L[line2][i+3:])
          name=word[0]
          occupancy = 1
          if name[0]=='*':
            name=name[1:]
            occupancy = 0
          if name[0]=='D':
            name='H'+name[1:]
          #print name
          anames[j]=name
          sx=[]
          for k in range(3):
            x=eval(word[k+1])
            while x<origin[k]:
              x=x+cell[k]
            while x>end[k]:
              x=x-cell[k]
            sx.append(str(x))
          #print name+'\t'+x+'\t'+y+'\t'+z
          acoords[j]=sx[0]+' '+sx[1]+' '+sx[2]
          #print acoords[j]
          j=j+1
          break
      i=i+1
    if flag1!=flag2: break
########################################################################
def readAtomicMovables(ifile,movables):
  ifile.seek(0)
  lines=ifile.readlines()
  
  searchterm='IONIC POSITIONS AND DISPLACEMENTS'
  
  j=0
  flag1=0
  flag2=0
  for line in lines: ## loop over lines of file 
    num_matches = string.find(line, searchterm)
    if num_matches>=0:
      flag1=1
    if flag1:
      if line[0]=='$' and line[1]=='$':
        flag2=1
        words=string.split(line[3:])
        name=words[0]
        movable = 1
        if name[0]=='*':
          movable = 0
        movables[j]=movable
        j=j+1
      else:
        if flag2:
          break      

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
    for i in range(0,natoms):
      print lineout[i]

########################################################################

searchterm1='IONIC POSITIONS AND FORCES'
searchterm2='Stepper Forces'

coords=[]
names=[]
movables=[]
for i in range(0,na):
  coords.append(' ')
  names.append(' ')
  movables.append(' ')

count_sets=0
for line in range(l): ## loop over lines of file 
  num_matches1 = string.find(L[line], searchterm1)
  num_matches2 = string.find(L[line], searchterm2)
  if num_matches1>=0 or num_matches2>=0 :
    modulus=count_sets%dump_freq_
    if( modulus==0 ):
      readAtomicPositions(line+1,line+na+2,names,coords)
      readAtomicMovables(input_,movables)
      if( dump_freq_<default_dump_freq_ ):
        filename=filename_+`count_sets`+".pdb"
        printPDBfile(na,names,coords,movables,filename)
    count_sets=count_sets+1

if( dump_freq_==default_dump_freq_ ):
  printPDBfile(na,names,coords,movables,"")
