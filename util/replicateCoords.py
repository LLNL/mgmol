# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#run:
#python replicateCoords.py coords.in mgmol_quench.cfg 2 2 2
import sys, string
coords=open(sys.argv[1],'r')
cfg=open(sys.argv[2],'r')
nx=eval(sys.argv[3])
ny=eval(sys.argv[4])
nz=eval(sys.argv[5])

nnn=nx*ny*nz

# list possible species in a dictionary with their count
species={'H':0,'He':0,'Li':0,'Be':0,'B':0,'C':0,'N':0,'O':0,'F':0,'Na':0,
         'Mg':0,'Al':0,'Si':0,'P':0,'S':0,'Cl':0,'K':0,'Ca':0,'D':0,'In':0,
         'Au':0,'Mo':0,'La':0}

lx=0.
ly=0.
lz=0.
for line in cfg:
  if "lx" in line:
    word=line.split('=')
    lx=eval(word[1])
  if "ly" in line:
    word=line.split('=')
    ly=eval(word[1])
  if "lz" in line:
    word=line.split('=')
    lz=eval(word[1])

count_atom = 0
for line in coords: ## loop over lines of file
  word=line.split()
  if len(word)>0:
    if len(word)>4:
      name1=word[0][0:1]
      name2=word[0][0:2]
      myspecies = 'none'
      if name1 in species.keys():
        myspecies = name1
      if name2 in species.keys():
        myspecies = name2
      if myspecies !='none':
        if len(word)>6:
          vx=eval(word[6])
          vy=eval(word[7])
          vz=eval(word[8])
        for i in range(nx):
          x=round(eval(word[2])+i*lx,12)
          for j in range(ny):
            y=round(eval(word[3])+j*ly,12)
            for k in range(nz):
              z=round(eval(word[4])+k*lz,12)
              count_atom = count_atom + 1
              name = myspecies + str(count_atom)
              sp=word[1]
              print (name.ljust(7),str(sp).rjust(3),str(x).rjust(16),str(y).rjust(16),str(z).rjust(16),end='')
              if len(word)>5:
                print (" ",word[5].rjust(3),end='')
              if len(word)>6:
                print (" ",str(vx).rjust(16),str(vy).rjust(16),str(vz).rjust(16),end='')
              print  ('\n',end='')
      else:
        print (line, end='')
    else:
      print (line,end='')
  else:
    print (line,end='')
