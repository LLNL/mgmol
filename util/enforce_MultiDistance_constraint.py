# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Python program to enforce constraint in mgmol input
#
# use: python enforce_MultiDistance_constraint.py coords.in alpha1 atom1a atom1b alpha2 atom2a atom2b ... distance
#-------------------------------------------------------------------------------
import sys, string
from math import sqrt
from numpy import *

##############################################
# count atoms/species
def getMassAtom(name):
  mass=0.
  if name[0:2] in spmass.keys():
     mass=spmass[name[0:2]]
  else:
     mass=spmass[name[0:1]]
  return mass

def getName(strn):
  if strn[0:2] in spmass.keys():
     return strn[0:2]
  else:
     return strn[0:1]

##############################################
# make one shake iteration
def enforce_constraint():
  sigma=(-1.)*distance_
  for k in range(nc_):
    x1=taup_[m_tau1_[k]]
    x2=taup_[m_tau2_[k]]
    dx = eval(x1[0])-eval(x2[0])
    dy = eval(x1[1])-eval(x2[1])
    dz = eval(x1[2])-eval(x2[2])
    d2 = dx*dx + dy*dy + dz*dz
    d = sqrt( d2 )
    sigma = sigma+eval(alpha_[k])*d

  print("#sigma=={}".format(sigma))
  if ( abs( sigma ) < tol_ ):
    return 1

  #make one shake iteration
  for p in range(na_):
  
    for e in range(3):
      a_[p][e]=0.
    
      for k in range(nc_):
        x1=taup_[m_tau1_[k]]
        x2=taup_[m_tau2_[k]]
        dx = eval(x1[0])-eval(x2[0])
        dy = eval(x1[1])-eval(x2[1])
        dz = eval(x1[2])-eval(x2[2])
        d  = sqrt( dx*dx + dy*dy + dz*dz )
        
        d1=0.
        d2=0.
        if u_name_[p]==name1_[k]:
          d1=1.
        if u_name_[p]==name2_[k]:
          d2=1.
        dd=d1-d2
        x1=tau_[m_tau1_[k]]
        x2=tau_[m_tau2_[k]]
        dtau=eval(x1[e])-eval(x2[e])
        a_[p,e]=a_[p,e]+eval(alpha_[k])*dd*dtau/d

  c=0.
  for p in range(na_):
    for e in range(3):
      c = c+a_[p,e] * a_[p,e] / u_mass_[p]

  vlambda = (-1.)*sigma/c
  for p in range(na_):
    x=taup_[u_taup_[p]]
    xx=zeros(3)
    for e in range(3):
      xx[e] = eval(x[e])+vlambda * a_[p,e] / u_mass_[p]
    taup_[u_taup_[p]]=[str(xx[0]),str(xx[1]),str(xx[2])]
  
  return 0

# list possible species in a dictionary with their mass
spmass={'H':1.0,'He':4.0,'Li':7.,'Be':9.,'B':10.8,'C':12.0,'N':14,'O':16.0,'F':19.0,'Na':23.,
        'Mg':24.3,'Al':27.0,'Si':28.1,'P':31.,'S':32.1,'Cl':35.5,'Ca':40.,'Ge':72.6,'D':2.0,
        'Zn':65.39}

input=open(sys.argv[1],'r')

#define constraint
i=2
tol_=1.e-8
alpha_=[]
name1_=[]
name2_=[]
while i<(len(sys.argv)-1):
  alpha_.append(sys.argv[i])
  name1_.append(sys.argv[i+1])
  name2_.append(sys.argv[i+2])
  print("#alpha={}, name1_={}, ,name2_={}".format(sys.argv[i],sys.argv[i+1],sys.argv[i+2]))
  i=i+3

distance_=eval(sys.argv[i])
print("#distance_={}".format(distance_))

nc_=len(alpha_)
print("#nc={}".format(nc_))

#list of unique names
u_name_=[]
u_mass_=[]

names_=[]

#read file and get data for constraint atoms
lines=input.readlines()

tau_=[]
names_=[]
masses=[]
found=0

for line in lines: ## loop over lines of file
  words=line.split()
  if len(words)>4:
    name = words[0]
    #print(name)
    if( words[0] in name1_ or words[0] in name2_ ):
      found=found+1
      names_.append(words[0])
      x=words[2]
      y=words[3]
      z=words[4]
      tau_.append([x,y,z])

na_=len(names_)
print("#na={}".format(na_))

#copy
taup_=tau_[:]

#setup

m_tau1_=[]
m_tau2_=[]

u_taup_=[]
for i in range(nc_):
  
  for ia in range(na_):
  
    if (name1_[i] == names_[ia] ):
      found = 0
      for j in range(len(u_name_)):
        if u_name_[j] == name1_[i] :
          found = 1
      if not found:
        u_name_.append(names_[ia])
        mass = getMassAtom(names_[ia])
        u_mass_.append(mass)
        u_taup_.append(ia)
      m_tau1_.append(ia)
    else:
      if (name2_[i] == names_[ia] ):
        found = 0
        for j in range(len(u_name_)):
          if u_name_[j] == name2_[i] :
            found = 1
        if not found:
          u_name_.append(names_[ia])
          mass = getMassAtom(names_[ia])
          u_mass_.append(mass)
          u_taup_.append(ia)
        m_tau2_.append(ia)

a_=zeros((na_,3))

##############################################

#print("#Old coordinates")
#for i in range(len(u_name_)):
#  print("#{}\t{}".format(u_name_[i],taup_[u_taup_[i]]))

conv=0
maxit=10
it=0
while conv==0:
  conv=enforce_constraint()
  it=it+1
  if it>maxit:
    break
    
#print("#New coordinates")
#for i in range(len(u_name_)):
#  print("#{}\t{}".format(u_name_[i],taup_[u_taup_[i]]))

if conv==1:
  #print("#write new file")
  for line in lines: ## loop over lines of file
    found=0
    words=line.split()
    if len(words)>3:
      for ia in range(len(names_)):
        if( words[0]==names_[ia] ):
          coord=taup_[ia]
          found=1
    if found==1:
      print(words[0].ljust(7),words[1].ljust(3), end ="")
      print("%10.6f  %10.6f  %10.6f"% (eval(coord[0]),eval(coord[1]),eval(coord[2])))
    else:
      print(line, end ="")

