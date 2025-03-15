#!/usr/bin/env python
import sys
import os
import subprocess
import string
import shutil

print("Test RestartEnergyAndForces...")

nargs=len(sys.argv)

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
for i in range(4,nargs-7):
  mpicmd = mpicmd + " "+sys.argv[i]
print("MPI run command: {}".format(mpicmd))

mgmol_exe = sys.argv[nargs-6]
test_exe = sys.argv[nargs-5]
input1 = sys.argv[nargs-4]
input2 = sys.argv[nargs-3]
coords = sys.argv[nargs-2]
print("coordinates file: %s"%coords)

#create links to potentials files
dst1 = 'pseudo.H_ONCV_PBE_SG15'
src1 = sys.argv[-1] + '/' + dst1

dst2 = 'pseudo.O_ONCV_PBE_SG15'
src2 = sys.argv[-1] + '/' + dst2

if not os.path.exists(dst1):
  print("Create link to %s"%dst1)
  os.symlink(src1, dst1)

if not os.path.exists(dst2):
  print("Create link to %s"%dst2)
  os.symlink(src2, dst2)

#run mgmol
command = "{} {} -c {} -i {}".format(mpicmd,mgmol_exe,input1,coords)
print("Run command: {}".format(command))

output = subprocess.check_output(command,shell=True)
lines=output.split(b'\n')

#analyse output
ref_energy=1.e18
for line in lines:
  if line.count(b'%%'):
    print(line)
    words=line.split()
    words=words[5].split(b',')[0]
    energy = words.decode()
  if line.count(b'achieved'):
    ref_energy=energy
    break

#run test
command = "{} {} -c {} -i {}".format(mpicmd,test_exe,input2,coords)
print("Run command: {}".format(command))
output = subprocess.check_output(command,shell=True)
lines=output.split(b'\n')

shutil.rmtree('WF')

test_energy=1.e18
l=-1
for line in lines:
  if line.count(b'Positions'):
    l=0
  if l>=0 and l<4:
    print(line)
    l=l+1
  if line.count(b'%%'):
    print(line)
    words=line.split()
    words=words[5].split(b',')[0]
    energy = words.decode()
  if line.count(b'Eks'):
    print(line)
    words=line.split()
    print(words)
    test_energy = words[2]
    break


tol = 1.e-6
diff=eval(test_energy)-eval(ref_energy)
print(diff)
if abs(diff)>tol:
  print("Energies differ: {} vs {} !!!".format(ref_energy,test_energy))
  sys.exit(1)

print("Test SUCCESSFUL!")
sys.exit(0)
