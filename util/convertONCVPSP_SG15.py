# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#convert qbox xml ONCV pseudopotentials into MGmol format
#usage example:
#python convertQboxONCVPSP.py O_ONCV_PBE-1.0.xml > pseudo.O_ONCV_PBE_SG15
import sys, string
import xml.etree.ElementTree as ET

filename=sys.argv[1]
tree = ET.parse(filename)
root = tree.getroot()
#for child in root:
#  print child.tag

description=root.find('description').text
pseudopot=root.find('norm_conserving_semilocal_pseudopotential')

symbol=root.find('symbol').text
atomic_number=eval(root.find('atomic_number').text)
mass=eval(root.find('mass').text)

valence_charge=eval(pseudopot.find('valence_charge').text)
mesh_spacing=pseudopot.find('mesh_spacing').text
h=eval(mesh_spacing)

local_pot=pseudopot.find('local_potential')
lpot=local_pot.text.split()

lmax=0
for projector in pseudopot.findall('projector'):
  l=int(projector.get('l'))
  if l>lmax:
    lmax=l

dij_factors=[]
for l in range(lmax+1):
  dij_factors.append([])

dijs=pseudopot.findall('d_ij')
for dij in dijs:
  l=eval(dij.get('l'))
  i=dij.get('i')
  j=dij.get('j')
  if i==j:
    v=dij.text
    dij_factors[l].append(v)

#print header
strings=description.split('\n')
for s in strings:
  print ('#',s)

print (filename.split('.')[0])
print ('#')
print ('color')
print ('#radii of balls and covalent bonds')
print ('-1.  -1.')
print ('# Nlcc flag')
print (0)
print ('# Atomic number')
print (atomic_number)
print ('# Atomic mass')
print (mass)
print ('# Number of valence electrons')
print (valence_charge)
print ('#Gaussian core charge parameter rc')
print (1.)
print ('# Number of potentials')
print (lmax+2)
print ('# l-value for state which is local, then type of potential format')
print (lmax+1, 3)
print ('# Local potential radius')
print (3.2)
print ('# Non-local potential radius')
print (3.2)

print ('# number of points in radial grid')
proj= pseudopot.find('local_potential')
ngpts=eval(proj.get('size'))
print (ngpts)

projs=[]
for l in range(lmax+1):
  projs.append([])

print ('# VANDERBILT-KLEINMAN-BYLANDER PROJECTORs')
print ('# l, nproj')
l=0
for dij in dij_factors:
  print (l,len(dij),end=" ")
  for d in dij:
    print (d,end=" ")
  print ('')
  l=l+1

for projector in pseudopot.findall('projector'):
  l=eval(projector.get('l'))
  i=projector.get('i')
  pot=projector.text
  s=pot.split()
  projs[l].append(s)

#print projectors
l=0
for proj in projs:
  print ('# l=',l)
  for i in range(ngpts):
    r=round(i*h,2)
    print (r,end=" ")
    for p in proj:
      print (eval(p[i]),end=" ")
    print ('')
  l=l+1
  
#print local potential
print ('# local')
for i in range(ngpts):
  r=round(i*h,2)
  print (r, lpot[i])
