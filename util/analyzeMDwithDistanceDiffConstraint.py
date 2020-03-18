# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#how to run in local directory:
#python analyzeMDwithDistanceDiffConstraint.py . O8203 H3041 N6678 > forces.dat
import sys, string, os
from math import sqrt, acos
from numpy import *
import matplotlib.pyplot as plt


kb_au=3.16678939e-06 # [Ha/K]
au2ps = 2.418885e-05

# list possible species in a dictionary with their mass
spmass={'H':1.0,'He':4.0,'Li':7.,'Be':9.,'B':10.8,'C':12.0,'N':14.,'O':16.0,'F':19.0,'Na':23.,
        'Mg':24.3,'Al':27.0,'Si':28.1,'P':31.,'S':32.1,'Cl':35.5,'Ca':40.,'Ge':72.6,'D':2.0}

avg_force=0.
running_avg_force=0.
corrected_avg_force=0.
avg_invsqrtz=0.

forces=[]
times=[]
running_forces=[]

nsteps=0
nravg=0
time=0.

def getMassAtom(name):
  one=name[0]
  mass=0.
  if one[0:2] in spmass.keys():
     mass=spmass[one[0:0+2]]
  else:
     mass=spmass[one[0]]
  return mass

def analyzeForces(filename,name0,name1,name2,dt):
  global avg_force
  global running_avg_force
  global corrected_avg_force
  global avg_invsqrtz
  global forces
  global nsteps
  global nravg
  global time
  
  m0=getMassAtom(name0)
  m1=getMassAtom(name1)
  m2=getMassAtom(name2)
  #print '#Masses:',m0,m1,m2

  coords=[]
  found=0
  found_force=0
  force=''

  r0=zeros(3)
  r1=zeros(3)
  r2=zeros(3)

  file=open(filename,'r')
  L1=file.readlines()
  for line in L1: ## loop over lines of file
    words=string.split(line)
    if len(words)>1:
      #find force on constraint
      if 'force' in words and 'constraint' in words and '=' in words:
        force=words[11]
        found_force=1
        #print 'force=',force
      if words[0][0:1]=='#' and found_force:
        name=words[1]
        if( name0==name or name1==name or name2==name):
          found=found+1
          x=words[2]
          y=words[3]
          z=words[4]
          coords.append(name+'\t'+x+'\t'+y+'\t'+z)
    if found==3 and found_force:
      for i in range(3):
        swords=string.split(coords[i])
        name=swords[0]
        x=eval(swords[1])
        y=eval(swords[2])
        z=eval(swords[3])
        if name==name0:
          r0=array([x,y,z])
        if name==name1:
          r1=array([x,y,z])
        if name==name2:
          r2=array([x,y,z])

      r01=r0-r1
      r12=r2-r1

      d01=sqrt((r01[0])**2+(r01[1])**2+(r01[2])**2)
      d12=sqrt((r12[0])**2+(r12[1])**2+(r12[2])**2)
      q=d01-d12
      #print 'Q=',q
      
      norm01=sqrt(inner(r01,r01))
      norm12=sqrt(inner(r12,r12))
      
      theta=acos(inner(r01,r12)/(norm01*norm12))
      z=1./m0+1./m2+4.*sin(0.5*theta)*sin(0.5*theta)/m1
      g=-1.*q*sin(theta)*sin(theta)/(z*z*m1*norm01*norm12)
      
      #print 'theta=',theta,', Z=',z,' G=',g
      corrected_force=eval(force)-kb_au*300.*g
      
      running_avg_force=(running_avg_force*nravg+eval(force))/(nravg+1)
      nravg=nravg+1

      times.append(time)
      forces.append(eval(force))
      
      print filename,'time=',time,' Q=',q,' force=',force,' corrected force=',corrected_force,' Z=',z,' running avg.= ',running_avg_force
      #print r01, r12
      time=time+dt*au2ps
      
      found=0
      coords=[]
      
      if time>1.:
        avg_force=avg_force+eval(force)
        
        running_forces.append(running_avg_force)

        invsqrtz=1./sqrt(z)
        avg_invsqrtz=avg_invsqrtz+invsqrtz
        corrected_avg_force=corrected_avg_force+corrected_force*invsqrtz
      
        nsteps=nsteps+1

  file.close()

#main
filesdir=sys.argv[1]
filenames=os.listdir(filesdir)

name0=sys.argv[2]
name1=sys.argv[3]
name2=sys.argv[4]
#print name0, name1, name2

inputs=[]
for filename in filenames:
  if 'md_run' in filename:
    inputs.append(filename)

inputs.sort()

inputs[0]
file=open(inputs[0],'r')
L1=file.readlines()
for line in L1: ## loop over lines of file
  word=string.split(line)
  if len(word)>1:
    if word[0]=='Timestep':
      dt=eval(word[5])

for filename in inputs:
  analyzeForces(filename,name0,name1,name2,dt)

nf=int(ceil(1./(dt*au2ps)))
ff=running_forces[-nf:]
maxF=max(float(v) for v in ff)
minF=min(float(v) for v in ff)
print '#maxF=',maxF
print '#minF=',minF

if nsteps>0:
  print '#Average force          =',avg_force/nsteps
  print '#Average corrected force=',corrected_avg_force/avg_invsqrtz
  print '#running_forces, spread last ',nf,' steps=',maxF-minF
  
skip=5
forcep = [ forces[i] for i in range(0, len(forces), skip)]
timesp = [ times[i] for i in range(0, len(forces), skip)]

skip=25
aves1 = [sum(forces[i-nf:i])/nf       for i in range(nf,   len(forces), skip)]
aves2 = [sum(forces[i-2*nf:i])/(2*nf) for i in range(2*nf, len(forces), skip)]
aves3 = [sum(forces[i-3*nf:i])/(3*nf) for i in range(3*nf, len(forces), skip)]

times1 = [ times[i] for i in range(nf, len(forces), skip)]
times2 = [ times[i] for i in range(2*nf, len(forces), skip)]
times3 = [ times[i] for i in range(3*nf, len(forces), skip)]

print '#time    force'
for i in range(len(forces)):
  print times[i],forces[i]

print '#Running average over 2 ps, last value: ',aves2[-1]
print '#Running average over 3 ps, last value: ',aves3[-1]


xmax=len(forces)*dt*au2ps
xmax=max(xmax,4.)
ymin=min(forces)
ymax=max(forces)

plt.figure(1)
plt.subplot(211)
plt.axis([0.,xmax,ymin,ymax])
plt.plot(timesp, forcep, 'go')
plt.ylabel('force (Ha/Bohr)')
#plt.show()

plt.subplot(212)
ymin=min(aves2)
ymax=max(aves2)
plt.axis([0.,xmax,ymin,ymax])
#plt.plot(times1, aves1, 'ro')
plt.plot(times2, aves2, 'ro')
plt.plot(times3, aves3, 'bo')
plt.ylabel('force (Ha/Bohr)')
plt.xlabel('time (ps)')
#plt.show()
plt.savefig('aves.png')

nn=nf/skip
max3=max(aves3[-1-nn:-1])
min3=min(aves3[-1-nn:-1])
print 'min aves3 over lat ps = ',min3
print 'max aves3 over lat ps = ',max3
print 'spread aves3 over lat ps = ',max3-min3
