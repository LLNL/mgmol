#!/usr/bin/env python
import sys
import os
import subprocess
import string

print("Test test_rho_restart...")

nargs=len(sys.argv)

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
for i in range(4,nargs-7):
  mpicmd = mpicmd + " "+sys.argv[i]
print("MPI run command: {}".format(mpicmd))

mgmol_exe = sys.argv[nargs-7]
test_exe = sys.argv[nargs-6]
input1 = sys.argv[nargs-5]
input2 = sys.argv[nargs-4]
input3 = sys.argv[nargs-3]
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

#run mgmol to generate initial ground state
command = "{} {} -c {} -i {}".format(mpicmd,mgmol_exe,input1,coords)
print("Run command: {}".format(command))

output = subprocess.check_output(command,shell=True)
lines=output.split(b'\n')

flag=0
for line in lines:
  if line.count(b'Run ended'):
    flag=1

if flag==0:
  print("Initial quench failed to complete!")
  sys.exit(1)

#run MD
command = "{} {} -c {}".format(mpicmd,mgmol_exe,input2)
print("Run command: {}".format(command))
output = subprocess.check_output(command,shell=True)
lines=output.split(b'\n')

flag=0
for line in lines:
  if line.count(b'Run ended'):
    flag=1

if flag==0:
  print("MD failed to complete!")
  sys.exit(1)

#run test
command = "{} {} -c {}".format(mpicmd,test_exe,input3)
print("Run command: {}".format(command))
output = subprocess.check_output(command,shell=True)
lines=output.split(b'\n')
for line in lines:
  print(line)

print("Test SUCCESSFUL!")
sys.exit(0)
