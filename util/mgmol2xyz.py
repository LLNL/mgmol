# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Python program to generate xyz files from mgmol output
#
# use:
# python mgmol2xyz.py mgmol_output > xyz
# python mgmol2xyz.py mgmol_output dump_freq
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
  if dump_freq_ == -1:
    dump_freq_ = default_dump_freq_
shiftx=0.
shifty=0.
shiftz=0.
if len(sys.argv)>3:
  if len(sys.argv)<6:
    print("Using shift require 6 arguments")
    exit(1)
  shiftx=eval(sys.argv[3])
  shifty=eval(sys.argv[4])
  shiftz=eval(sys.argv[5])

if len(sys.argv)>3:
  jobname=sys.argv[3]
else:
  jobname=''

lines=input_.readlines()
l=len(lines)  ## no of lines in file

# list possible species in a dictionary with their count
species={'H':0,'D':0,'He':0,'Li':0,'Be':0,'B':0,'C':0,'N':0,'O':0,'F':0,'Na':0,
         'Mg':0,'Al':0,'Si':0,'P':0,'S':0,'Cl':0,'Ca':0,'In':0,'Au':0,'X':0}

ang2bohr=1.8897269
bohr2ang=1./ang2bohr

na=outputTools.countNumAtoms(input_)
nmlwc=outputTools.countNumMLWC(input_)
#print("N atoms: {}".format(na))
#print("Number of MLWC: {}".format(nmlwc))

searchterm='Origin'
ox=0.
oy=0.
oz=0.
for line in lines:
  if line.count(searchterm):
    words=line.split()
    ox=eval(words[2])
    oy=eval(words[3])
    oz=eval(words[4])
origin=array([ox+shiftx,oy+shifty,oz+shiftz])

searchterm='Dimension'
lx=0.
ly=0.
lz=0.
for line in lines:
  if line.count(searchterm):
    words=line.split()
    lx=eval(words[2])
    ly=eval(words[3])
    lz=eval(words[4])

cell=array([lx,ly,lz])
end=array([lx+shiftx,ly+shifty,lz+shiftz])

########################################################################
def readAtomicPositions(first_line,last_line,anames,acoords):
  j=0
  flag1=0
  flag2=0
  #print("loop starting at {}".format(first_line+1))
  for line2 in range(first_line,last_line):
    i=0
    for c in lines[line2]:
      flag1=0
      if c=='#': 
        if lines[line2][i+1]=='#':
          flag1=1
          flag2=1
          word=lines[line2][i+3:].split()
          name=word[0]
          occupancy = 1
          if name[0]=='*':
            name=name[1:]
            occupancy = 0
          if name[0]=='D':
            name='H'+name[1:]
          x=word[1]
          y=word[2]
          z=word[3]
          #print(name+'\t'+x+'\t'+y+'\t'+z)
          anames[j]=name
          acoords[j]=x+' '+y+' '+z
          j=j+1
          break
      i=i+1
    if flag1!=flag2: break

########################################################################
def readMLWC(first_line,last_line,anames,acoords):
  count=0
  flag1=0
  flag2=0
  #print("loop starting at {}".format(first_line+1))
  for line2 in range(first_line,last_line):
    i=0
    for c in lines[line2]:
      flag1=0
      if c=='&': 
        if lines[line2][i+1]=='&' and 'Anderson' not in lines[line2]:
          flag1=1
          flag2=1
          word=lines[line2].split()
          x=word[2]
          y=word[3]
          z=word[4]
          #print(x+'\t'+y+'\t'+z)
          anames[count]='X'+str(count)
          acoords[count]=x+' '+y+' '+z
          count=count+1
          break
      i=i+1
    if flag1!=flag2: break

########################################################################
def makeXYZobject(natoms,anames,acoords,xyzobject):
  
  #atom name
  for i in range(natoms):
    if anames[i][0:2] in species.keys():
      name=anames[i][0:2]
    else:
      name=anames[i][0:1]
    name=name.ljust(4)
    xyzobject.append(name)

  #x,y,z
  for i in range(natoms): 
    word=acoords[i].split()
    xyzobject[i]=xyzobject[i]+'    '
    for j in range(3):
      x=eval(word[j])
      while x<origin[j]:
        x=x+cell[j]
      while x>end[j]:
        x=x-cell[j]
      x = "%.3f" % (x*bohr2ang)
      x=x.rjust(8)
      xyzobject[i]=xyzobject[i]+x

########################################################################
def printXYZfile(natoms,anames,acoords,filename,jobname=''):
  lineout=[]
  makeXYZobject(natoms,anames,acoords,lineout)

  if filename:
    output = open(filename,'w')
    output.write(str(natoms))
    output.write('\n')

    if jobname:
        output.write(jobname)
        output.write('\n')

    for i in range(0,natoms):
      output.write(lineout[i])
      output.write('\n')
  else:
    print(natoms)
    print("")
    for i in range(0,natoms):
      print(lineout[i])

########################################################################

searchterm1='IONIC POSITIONS AND FORCES'
searchterm2='Stepper Forces'
searchterm3='Orbitals centers and spreads'

acoords=[]
anames=[]
for i in range(0,na):
  acoords.append(' ')
  anames.append(' ')

wcoords=[]
wnames=[]
for i in range(0,nmlwc):
  wcoords.append(' ')
  wnames.append(' ')

#count number of frames
nframes=0
for line in lines:
  if( line.find('IONIC CONFIGURATION ENERGY')>0 ):
    nframes=nframes+1
#print("Number of frames: {}".format(nframes))

count_sets=0
for line in range(l): ## loop over lines of file 
  num_matches1 = lines[line].find(searchterm1)
  num_matches2 = lines[line].find(searchterm2)
  num_matches3 = lines[line].find(searchterm3)

  if num_matches3>=0 and nmlwc>0:
    readMLWC(line+1,line+nmlwc+2,wnames,wcoords)

  if num_matches1>=0 or num_matches2>=0 :
    modulus=count_sets%dump_freq_
    #print("Count = {}".format(count_sets))
    if( modulus==0 or count_sets==(nframes-1) ):
      #print("Read frame {}".format(count_sets))
      readAtomicPositions(line+1,line+na+2,anames,acoords)
      if( dump_freq_<default_dump_freq_ ):
        filename=filename_+str(count_sets)+".xyz"
        printXYZfile(na+nmlwc,anames+wnames,acoords+wcoords,filename, jobname)
    count_sets=count_sets+1
  

names =anames +wnames
coords=acoords+wcoords
  
if( dump_freq_==default_dump_freq_ ):
  printXYZfile(na+nmlwc,names,coords,"")
