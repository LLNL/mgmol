#!/usr/bin/env python
import sys
import os
import subprocess
import string

print("Test LBFGS...")

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
dst = 'pseudo.Si_ONCV_PBE_SG15'
src = sys.argv[-1] + '/' + dst

if not os.path.exists(dst):
  print("Create link to %s"%dst)
  os.symlink(src, dst)

#run MGmol
command = "{} {} -c {} -i {} -l {}".format(mpicmd,exe,inp1,coords,lrs)
print(command)
output = subprocess.check_output(command,shell=True)

#analyse mgmol standard output
lines=output.split(b'\n')

for line in lines:
  if line.count(b'CONFIG'):
    print(line)
  if line.count(b'lbfgs'):
    print(line)
  if line.count(b'max') and line.count(b'movable'):
    print(line)

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

command = "{} {} -c {} -l {}".format(mpicmd,exe,inp2,lrs)
print(command)
output = subprocess.check_output(command,shell=True)

#analyse mgmol standard output
lines=output.split(b'\n')

print("Check LBFGS convergence...")
tol = 4.e-4
force = 1.e6
for line in lines:
  if line.count(b'%%') and line.count(b'ION'):
    print(line)

  if line.count(b'max') and line.count(b'movable'):
    print(line)
    words=line.split()
    
    force=eval(words[6])

print("Last max force = {}".format(force))
os.remove('wave.out')

if force>tol:
  print("Force larger than tol {}".format(tol))
  sys.exit(1)

sys.exit(0)
