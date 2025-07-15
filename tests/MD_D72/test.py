#!/usr/bin/env python
import sys
import os
import subprocess
import string
import shutil

print("Test MD...")

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
command = "{} {} -c {} -i {} -l {}".format(mpicmd,exe,inp1,coords,lrs)
print("Run command: {}".format(command))
output1 = subprocess.check_output(command,shell=True)
lines=output1.split(b'\n')

#analyse output of quench
for line in lines:
  num_matches = line.count(b'%%')
  if num_matches:
    print(line)

#run MD
command = "ls -ld snapshot* | awk '{ print $9 }' | tail -n1"
print(command)
restart_file = subprocess.check_output(command,shell=True)
restart_file=str(restart_file[:-1],'utf-8')
print(restart_file)

try:
  os.symlink(restart_file, 'wave.out')
except FileExistsError:
  os.remove('wave.out')
  os.symlink(restart_file, 'wave.out')

command = "{} {} -c {}".format(mpicmd,exe,inp2)
output2 = subprocess.check_output(command,shell=True)

#remove created files
shutil.rmtree(restart_file)
os.remove('wave.out')

#analyse mgmol standard output
lines=output2.split(b'\n')

print("Check energy conservation...")
tol = 1.e-2
energy = 0.
count = 0
for line in lines:
  num_matches1 = line.count(b'Total')
  num_matches2 = line.count(b'Energy')
  if num_matches1 and num_matches2:
    print(line)
    count=count+1
    words=line.split()
    
    energy=eval(words[2])
    if count==1:
      first_energy=energy

    if count>1 and abs(energy-first_energy)>tol:
      print("ERROR Energy = {} != {}".format(energy,first_energy))
      sys.exit(1)

sys.exit(0)
