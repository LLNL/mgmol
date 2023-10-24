# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# How to run in a python script:
#   sys.path.append('/home/q8j/Documents/q8j/GIT/MGmol/util')
#   import analyzeMDwithDistanceDiffConstraint as analyze
#   ...
#   analyze.runAnalyzeForces(idir,name0,name1,name2)

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
  mass=0.
  if name[0:2] in spmass.keys():
     mass=spmass[name[0:0+2]]
  else:
     mass=spmass[name[0]]
  return mass

def analyzeForces(filename,name0,name1,name2,dt):
  global avg_force
  global running_avg_force
  global corrected_avg_force
  global avg_invsqrtz
  global forces
  global times
  global running_forces
  global nsteps
  global nravg
  global time

  print("analyzeForces with arguments {}, {}, {}, {}".format(filename, name0,name1,name2,dt))

  m0=getMassAtom(name0)
  m1=getMassAtom(name1)
  m2=getMassAtom(name2)
  print("#Masses: {}, {}, {}".format(m0,m1,m2))

  coords=[]
  found=0
  found_force=0
  force=''

  r0=zeros(3)
  r1=zeros(3)
  r2=zeros(3)

  f=open(filename,'r')
  lines=f.readlines()
  for line in lines: ## loop over lines of file
    words=line.split()
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
        swords=coords[i].split()
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

      tol = 1.e-4
      delta = corrected_force-eval(force)
      if abs(delta)>tol:
        print("{}, time={}, Q={}, force={}, corrected force={}, Z={}, running avg.= {}".format( \
              filename,time,q,force,corrected_force,z,running_avg_force))
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

  f.close()

###############################################################################
# main script
###############################################################################


def runAnalyzeForces(filesdir,name0,name1,name2):

  global avg_force
  global running_avg_force
  global corrected_avg_force
  global avg_invsqrtz
  global forces
  global times
  global running_forces
  global nsteps
  global nravg
  global time

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


  filenames=os.listdir(filesdir)

  #screen out files that are not mgmol outputs
  inputs=[]
  for filename in filenames:
    if '.out' in filename:
      inputs.append(filesdir+"/"+filename)

  inputs.sort()

  #extract dt from first output file in MD run
  inputs[0]
  f=open(inputs[0],'r')
  lines=f.readlines()
  for line in lines: ## loop over lines of file
    word=line.split()
    if len(word)>1:
      if word[0]=='Timestep':
        dt=eval(word[5])

  #analyze all outputs
  for filename in inputs:
    analyzeForces(filename,name0,name1,name2,dt)

  #number of steps/ps
  nf=int(ceil(1./(dt*au2ps)))
  print("nf = {}".format(nf))

  ff=running_forces[-nf:]
  maxF=max(float(v) for v in ff)
  minF=min(float(v) for v in ff)
  print("#maxF={}".format(maxF))
  print("#minF={}".format(minF))

  if nsteps>0:
    print("#Average force          ={}".format(avg_force/nsteps))
    print("#Average corrected force={}".format(corrected_avg_force/avg_invsqrtz))
    print("#running_forces, spread last {}, steps={}".format(nf,maxF-minF))
  
  skip=5
  forcep = [ forces[i] for i in range(0, len(forces), skip)]
  timesp = [ times[i] for i in range(0, len(forces), skip)]
  plt.figure(1)
  xmax=len(forces)*dt*au2ps
  ymin=min(forces)
  ymax=max(forces)
  plt.axis([0.,xmax,ymin,ymax])

  fig=plt.figure(1)
  plt.subplot(211)
  plt.plot(timesp, forcep, 'go')
  plt.ylabel('force (Ha/Bohr)')
  colors=['ro','bo','go','yo','mo']
  icolor=0
  #compute running averages
  skip=25
  for interval in range(1,6):
    aves = [sum(forces[i-interval*nf:i])/(interval*nf) for i in range(interval*nf,len(forces), skip)]

    timesp = [ times[i] for i in range(interval*nf, len(forces), skip)]
    forcep = [ forces[i] for i in range(interval*nf, len(forces), skip)]

    print("#Average over last {} ps: {}".format(interval,aves[-1]))

    sigma = [sum((forces[i-interval*nf:i]-aves[-1])*(forces[i-interval*nf:i]-aves[-1]))/(interval*nf-1) for i in range(interval*nf, len(forces), skip)]

    print("#Sigma over last {} ps: {}".format(interval,math.sqrt(sigma[-1])))

    plt.subplot(212)
    plt.plot(timesp, aves, colors[icolor])
    icolor=icolor+1
    plt.ylabel('running average (Ha/Bohr)')
    plt.xlabel('time (ps)')
    plt.axis([0.,xmax,ymin,ymax])

  plt.savefig(filesdir+'_aves.png')
  plt.close(fig)

  nn=int(nf/skip)
  print(nn)
  max3=max(aves[-1-nn:-1])
  min3=min(aves[-1-nn:-1])
  print("min aves over last ps = {}".format(min3))
  print("max aves over last ps = {}".format(max3))
  print("min-max spread aves3 over last ps = {}".format(max3-min3))

###############################################################################
# main script
###############################################################################

#target directory containing multiple segments of a single run
#filesdir=sys.argv[1]
#filenames=os.listdir(filesdir)

#read names of 3 atoms involved in constraint
#name0=sys.argv[2]
#name1=sys.argv[3]
#name2=sys.argv[4]

#runAnalyzeForces(filesdir,name0,name1,name2)
