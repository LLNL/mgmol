# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#run:
#python replicateInput.py filename.in nx ny nz > new_filename.in
import sys, string
input=open(sys.argv[1],'r')
nx=eval(sys.argv[2])
ny=eval(sys.argv[3])
nz=eval(sys.argv[4])

# list possible species in a dictionary with their count
species={'H':0,'He':0,'Li':0,'Be':0,'B':0,'C':0,'N':0,'O':0,'F':0,'Na':0,
         'Mg':0,'Al':0,'Si':0,'P':0,'S':0,'Cl':0,'K':0,'Ca':0,'D':0,'In':0}


count=0
count_atom = 0
reached_atoms = 0
nsp = 0
Li=input.readlines()
l=len(Li)  ## no of lines in file
for line in range(l): ## loop over lines of file
  word=string.split(Li[line])
  if len(word)>0:
    if word[0][0:1]!='#':
      count=count+1
      if count==3:
        nnx=nx*eval(word[0])
        nny=ny*eval(word[1])
        nnz=nz*eval(word[2])
        print nnx,'\t',nny,'\t',nnz
      elif  count==2:
        b1=eval(word[0])
        b2=eval(word[1])
        b3=eval(word[2])
        e1=eval(word[3])
        e2=eval(word[4])
        e3=eval(word[5])
        d1=e1-b1
        d2=e2-b2
        d3=e3-b3
        print b1,'\t',b2,'\t',b3,'\t',
        print e1+(nx-1)*d1,'\t',e2+(ny-1)*d2,'\t',e3+(nz-1)*d3
      elif  count==4:
        nsp=eval(word[0])
        print nsp
      elif  count==32+nsp:
        ne=eval(word[0])
        print ne*nx*ny*nz
      elif  count==33+nsp:
        ne=eval(word[0])
        print ne*nx*ny*nz,'  ',eval(word[1])
      elif  count==34+nsp:
        ni=eval(word[0])
        print ni*nx*ny*nz
      else:
        if count>40+nsp and len(word)>4:
          name1=word[0][0:1]
          name2=word[0][0:2]
          myspecies = 'none'
          if name1 in species.keys():
            myspecies = name1
          if name2 in species.keys():
            myspecies = name2
          if myspecies !='none':
            reached_atoms = 1
            for i in range(nx):
              x=eval(word[2])+i*d1
              for j in range(ny):
                y=eval(word[3])+j*d2
                for k in range(nz):
                  z=eval(word[4])+k*d3
                  count_atom = count_atom + 1
                  name = myspecies + str(count_atom)
                  sp=word[1]
                  print name.ljust(7),str(sp).rjust(3),str(x).rjust(16),str(y).rjust(16),str(z).rjust(16),
                  if len(word)>5: print word[5].rjust(3)
                  else: print  '\n',
          else:
            print Li[line],  
        else:
          if ( len(word)==4 or len(word)==3 ) and reached_atoms > 0:
            if len(word)==4:
              for i in range(nx):
                x=eval(word[0])+i*d1
                for j in range(ny):
                  y=eval(word[1])+j*d2
                  for k in range(nz):
                    z=eval(word[2])+k*d3
                    print str(x).ljust(16),str(y).ljust(16),str(z).ljust(16),word[3].rjust(4)
            else:
              for i in range(nx):
                x=eval(word[0])+i*d1
                for j in range(ny):
                  y=eval(word[1])+j*d2
                  for k in range(nz):
                    z=eval(word[2])+k*d3
                    print str(x).ljust(16),str(y).ljust(16),str(z).ljust(16)
          else:
            print Li[line],
    else:
      print Li[line],
  else:
    print Li[line],
