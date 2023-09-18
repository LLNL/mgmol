# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#compute domain boundaries to include atoms in PDB file with a buffer region
#output in Bohr
import sys, string
myfile=open(sys.argv[1],'r')
lines=myfile.readlines()

bufsize=eval(sys.argv[2])

ang2bohr=1.8897269

minx=1000.
miny=1000.
minz=1000.
maxx=-1000.
maxy=-1000.
maxz=-1000.
for line in lines:
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

#print lower left corner
print ('ll=(',ang2bohr*minx-bufsize,ang2bohr*miny-bufsize,ang2bohr*minz-bufsize,')')

#print domain dimensions
print ('lenth=(',ang2bohr*(maxx-minx)+2.*bufsize, \
                ang2bohr*(maxy-miny)+2.*bufsize, \
                ang2bohr*(maxz-minz)+2.*bufsize,')')
