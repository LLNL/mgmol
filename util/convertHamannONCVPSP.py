# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#convert ONCV pseudopotentials from Hamann's code into MGmol format
#usage example:
#python convertHamannONCVPSP.py 08_O.out > pseudo.O_ONCV_PBE
import sys, string

filename=sys.argv[1]
myfile=open(sys.argv[1],'r')
lines=myfile.readlines()

# list possible species in a dictionary with their mass
spmass={'H':1.0,'He':4.0,'Li':7.,'Be':9.,'B':10.8,'C':12.0,'N':14.,'O':16.0,
        'F':19.0,'Na':23.,'Mg':24.3,'Al':27.0,'Si':28.1,'P':31.,'S':32.1,
        'Cl':35.5,'Ca':40.,'Ge':72.6,'D':2.0}

#general info
print("#Generated by Hamann's code")

active=False
count=0
for line in lines:
  if count<2:
    print("#"+line,end='')
  words=line.split()
  if len(words)>0:
    if words[0]=='#':
      active=True
    if words[0]=='Reference':
      break
    else:
      if active:
        print("#"+line,end='')
  count=count+1

print("#Reference:")
print("#D. R. Hamann, Phys. Rev. B 88, 085117 (2013)")

#header paramaters
dij_factors=[]

active = False
step=0
mesh_spacing=0.01
nprojs=[]
ip=0
il=0
lloc=-1
for line in lines:
  words=line.split()
  if len(words)>1:
    if words[1]=='PSPCODE8':
      active = True
  if active:
    if step==1:
      print("# Short description of the species type. One line only!")
      description=words[0]+' '+words[1]
      symbol=words[0]
      print(description)
      print("#color")
      print("Undefined")
      print("#radii of balls and covalent bonds")
      print("0. 0.")
      print("# Nlcc flag")
      print("0")
    if step==2:
      print("# Atomic number")
      print(words[0])
      print("# Atomic mass")
      print(spmass[symbol])
      print("# Number of valence electrons")
      print(words[1])
      print("# Gaussian core charge parameter rc")
      print("1.")
    if step==3:
      lmax=eval(words[2])
      lloc=eval(words[3])
      for l in range(lmax+1):
        dij_factors.append([])
      print("# Number of potentials")
      print (lmax+2)
      print ('# l-value for state which is local, then type of potential format')
      print (lmax+1, 2)
      print ("# Local potential radius")
      print (3.2)
      print ("# Non-local potential radius")
      print (3.2)
      print ("# number of points in radial grid")
      n=eval(words[4])
      print (n)
    if step==5:
      print ("# VANDERBILT-KLEINMAN-BYLANDER PROJECTORS")
      print ("# l, nproj")
      for l in range(lmax+1):
        nprojs.append(eval(words[l]))
      for l in range(lmax+1,lloc+1):
        nprojs.append(0)
      nprojs[lloc]=1

    if step==7: #coefficients
      print(il,' ',nprojs[il],end='  ')
      for l in range(nprojs[il]):
        coeff=float(words[l+1].replace("D","E"))
        print(" {}".format(coeff),end='')
      print()
      ip=0 
 
    if step==8: #read potentials
      #print(line)
      ip = ip+1
      if ip==n:
        ip=0
        step=7
        if il==lmax:
          break
        il=il+1
    else:
      step = step + 1

#now read potential functions
active=False
count=0
iread=-1
for line in lines:
  words=line.split()
  if len(words)>1:
    if words[1]=='PSPCODE8':
      active=True
  if active:
    #print (line)
    count=count+1
    if count%(n+1)==8:
      l=eval(words[0])
      print("# radial grid, and projecrtors for l={}".format(l))
    else:
      if count>7:
        r=float(words[1].replace("D","E"))
        print("{0:3.2f}".format(r),end='')
        vals=[]
        for i in range(nprojs[l]):
          val=float(words[i+2].replace("D","E"))
          print("  {}".format(val),end='')
        print('\n',end='')

sys.exit()

