#!/usr/bin/env python
import sys
import os
import subprocess
import string

print("Test continuum solvent...")

nargs=len(sys.argv)

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
for i in range(4,nargs-5):
  mpicmd = mpicmd + " "+sys.argv[i]
print("MPI run command: {}".format(mpicmd))

exe = sys.argv[nargs-5]
inp = sys.argv[nargs-4]
coords = sys.argv[nargs-3]
print("coordinates file: %s"%coords)
lrs = sys.argv[-2]

#create links to potentials files
dst1 = 'pseudo.O_ONCV_PBE_SG15'
dst2 = 'pseudo.H_ONCV_PBE_SG15'
src1 = sys.argv[nargs-1] + '/' + dst1
src2 = sys.argv[nargs-1] + '/' + dst2

if not os.path.exists(dst1):
  print("Create link to %s"%dst1)
  os.symlink(src1, dst1)
if not os.path.exists(dst2):
  print("Create link to %s"%dst2)
  os.symlink(src2, dst2)

#run mgmol
command = "{} {} -c {} -i {}".format(mpicmd,exe,inp,coords)
print("Run command: {}".format(command))
output = subprocess.check_output(command,shell=True)

os.remove(dst1)
os.remove(dst2)

#analyse mgmol standard output
#make sure energy is correct
lines=output.split(b'\n')

energy = 0.
for line in lines:
  if line.count(b'%%'):
    print(line)
    words = line.split()
    energy = eval(words[5][:-1])

target_energy = -17.1826
print("Last energy = {}".format(energy))
if abs(target_energy-energy) > 1.e-3:
  print("Incorrect energy")
  sys.exit(1)

sys.exit(0)
