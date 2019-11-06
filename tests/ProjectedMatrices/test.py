#!/usr/bin/env python
import sys
import os
import subprocess
import string

def getData(lines, energies, meanf):
  for line in lines:
      if line.count(b'%%'):
        print(line)
        words=line.split()
        energy = words[5].decode()
        energies.append(energy.replace(',', ''))
      if line.count(b'mean') and line.count(b'movable'):
        words=line.decode("ascii").split(",")
        meanf.append(words[0].split("(")[-1])
        meanf.append(words[1])
        meanf.append(words[2].split(")")[0])
        print("{} {} {}".format(meanf[0],meanf[1],meanf[2]))
        break

print("Test ProjectedMatrices...")

nargs=len(sys.argv)

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
for i in range(4,nargs-6):
  mpicmd = mpicmd + " "+sys.argv[i]
print("MPI run command: {}".format(mpicmd))

exe = sys.argv[nargs-6]
inp1 = sys.argv[nargs-5]
inp2 = sys.argv[nargs-4]
coords = sys.argv[nargs-3]
print("coordinates file: %s"%coords)
lrs = sys.argv[-2]

#create links to potentials files
dst = 'pseudo.D_ONCV_PBE_SG15'
src = sys.argv[-1] + '/' + dst

if not os.path.exists(dst):
  print("Create link to %s"%dst)
  os.symlink(src, dst)

#run quench
command1 = "{} {} -c {} -i {} -l {}".format(mpicmd,exe,inp1,coords,lrs)
print("Run command: {}".format(command1))
output1 = subprocess.check_output(command1,shell=True)
lines1=output1.split(b'\n')

#analyse output of quench1
energies1=[]
meanf1=[]
getData(lines1, energies1, meanf1)

#run 2nd quench
command2 = "{} {} -c {} -i {} -l {}".format(mpicmd,exe,inp2,coords,lrs)
print("Run command: {}".format(command2))
output2 = subprocess.check_output(command2,shell=True)
lines2=output2.split(b'\n')

#analyse output of quench2
energies2=[]
meanf2=[]
getData(lines2, energies2, meanf2)

if len(energies1) != len(energies2):
  print("Energies1:")
  print(energies1)
  print("Energies2:")
  print(energies2)
  sys.exit(1)

print("Check energies...")
tol = 1.e-6
for i in range(len(energies1)):
  energy1=eval(energies1[i])
  energy2=eval(energies2[i])
  diff=energy2-energy1
  if abs(diff)>tol:
    print("Energies differ: {} vs {} !!!".format(energies1[i],energies2[i]))
    sys.exit(1)

print("Check mean forces...")
if (eval(meanf2[0])-eval(meanf1[0]))>tol:
  print("mean F values in x-direction: {}, {}".format(meanf1[0],meanf2[0]))
  sys.exit(1)
if (eval(meanf2[1])-eval(meanf1[1]))>tol:
  print("mean F values in y-direction: {}, {}".format(meanf1[1],meanf2[1]))
  sys.exit(1)
if (eval(meanf2[2])-eval(meanf1[2]))>tol:
  print("mean F values in z-direction: {}, {}".format(meanf1[2],meanf2[2]))
  sys.exit(1)

print("Test SUCCESSFUL!")
sys.exit(0)
