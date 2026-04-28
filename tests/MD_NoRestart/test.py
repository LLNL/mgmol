#!/usr/bin/env python
import sys
import os
import subprocess
import string
import shutil

print("Test MD without restart...")

nargs=len(sys.argv)

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
for i in range(4,nargs-4):
  mpicmd = mpicmd + " "+sys.argv[i]
print("MPI run command: {}".format(mpicmd))

exe = sys.argv[nargs-4]
inp = sys.argv[nargs-3]
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

#run MD
command = "{} {} -c {} -i {}".format(mpicmd,exe,inp,coords)
output = subprocess.check_output(command,shell=True)

#analyse mgmol standard output
lines=output.split(b'\n')

print("Check energy conservation...")
tol = 1.e-2
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

for line in lines:
  if line.count(b'Kinetic'):
    print(line)
    words=line.split()
    kinetic = eval(words[2])
    if kinetic<0.0004:
      print("ERROR Expect larger kinetic energy!")
      sys.exit(1)

if count < 4:
  print("Expect to finfd 4 energies!")
  sys.exit(1)

sys.exit(0)
