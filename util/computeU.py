# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
import numpy as np
from scipy.linalg import sqrtm

def format(filename):

    with open(filename, 'r') as f:

        lines = f.readlines()

    data = np.zeros(len(lines))

    for i in range(len(lines)):

        data[i] = float(lines[i].split('=')[1])
            
    numst = int(np.sqrt(len(data)))

    data = data.reshape(numst, numst)

    return data

if __name__ == '__main__':

    phi = np.loadtxt('phi.out')
    psi = np.loadtxt('psi.out')
    
    phi = np.flipud(phi)
    psi = np.flipud(psi)
    
    O = np.dot(psi, phi.transpose())
    O = O / O[0][0]
    
    OOT = np.dot(O, O.transpose())
    
    fac = np.linalg.inv(sqrtm(OOT))
    
    U = np.dot(fac, O)
    
    Uc = format('U.out')
    Oc = format('O.out')
    OOTc = format('OOT.out')
    Dc = format('D.out')
    Vc = format('V.out')
    VDc = format('VD.out')
    VDVTc = format('VDVT.out')
    Ec = np.loadtxt('eigenvalues.out')
    
    #print O
    
    #print oot
    
    #print fac
    
    print U
    print Uc
    print U - Uc
    print O - Oc
    
    import pdb
    pdb.set_trace()
