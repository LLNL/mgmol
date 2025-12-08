#!/usr/bin/env python
import sys
import os
import subprocess
import string

print("Test SpreadPenalty...")

nargs=len(sys.argv)

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
for i in range(4,nargs-5):
  mpicmd = mpicmd + " "+sys.argv[i]
print("MPI run command: {}".format(mpicmd))

exe = sys.argv[nargs-5]
inp = sys.argv[nargs-4]
coords = sys.argv[nargs-3]
lrs = sys.argv[nargs-2]
print("coordinates file: %s"%coords)

#create links to potentials files
dst = 'pseudo.O_ONCV_PBE_SG15'
src = sys.argv[-1] + '/' + dst
if not os.path.exists(dst):
  print("Create link to %s"%dst)
  os.symlink(src, dst)

dst = 'pseudo.H_ONCV_PBE_SG15'
src = sys.argv[-1] + '/' + dst
if not os.path.exists(dst):
  print("Create link to %s"%dst)
  os.symlink(src, dst)

#run MGmol
command = "{} {} -c {} -i {} -l {}".format(mpicmd,exe,inp,coords,lrs)
output = subprocess.check_output(command,shell=True)

#analyse mgmol standard output
lines=output.split(b'\n')

#print number MPI tasks and threads used
for i in range(50):
  line = lines[i]
  if line.count(b'active') and line.count(b'thread'):
    print(line)
  if line.count(b'MPI') and line.count(b'tasks'):
    print(line)

print("Check energy convergence and spread values...")

check_spread_init = True
found = False
spread = -1.
energy = 0.
for line in lines:
  if line.count(b'spread') and line.count(b'='):
    words=line.split()
    spread = eval(words[3])
    print("spread = {}".format(spread))
    if check_spread_init:
      if spread<3.:
        print("Test failed: first spread too small")
        sys.exit(1)

      check_spread_init = False
  if line.count(b'%%'):
    words=line.split(b',')
    words=words[0].split()
    energy = eval(words[-1])
    print("energy = {}".format(energy))

#make sure target spread is reached within 1%
spread_target = 1.5
tol = spread_target*0.01
if (spread-spread_target) > tol:
  print("Test failed: last spread too large!")
  sys.exit(1)

#we tolerate an energy difference since the initial wave functions
#are very delocalized and the spread penalty remains active all along
energy_ref = -17.1660
tol = 5.e-4
if abs(energy-energy_ref) > tol:
  print("Test failed: last energy value incorrect!")
  sys.exit(1)

sys.exit(0)
