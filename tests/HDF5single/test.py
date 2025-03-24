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

mgmol_exe = sys.argv[nargs-5]
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

#run mgmol to generate initial ground state
command = "{} {} -c {} -i {}".format(mpicmd,mgmol_exe,input1,coords)
print("Run command: {}".format(command))

output = subprocess.check_output(command,shell=True)
lines=output.split(b'\n')

#run MD
command = "{} {} -c {} -i {}".format(mpicmd,mgmol_exe,input2,coords)
print("Run command: {}".format(command))
output = subprocess.check_output(command,shell=True)
lines=output.split(b'\n')

os.remove('wf.h5')

print("Check energy conservation...")
tol = 1.e-4
energy = 0.
count = 0
for line in lines:
  if line.count(b'Total') and line.count(b'Energy'):
    print(line)
    count=count+1
    words=line.split()

    energy=eval(words[2])
    if count==1:
      first_energy=energy

    if count>1 and abs(energy-first_energy)>tol:
      print("ERROR Energy = {} != {}".format(energy,first_energy))
      sys.exit(1)

if count<4:
  print("ERROR needs 4 energy values for checking conservation!")
  sys.exit(1)

os.remove('wf_md.h5')

print("Test SUCCESSFUL!")
sys.exit(0)
