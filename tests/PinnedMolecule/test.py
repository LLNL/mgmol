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

exe = sys.argv[nargs-5]
inp1 = sys.argv[nargs-4]
inp2 = sys.argv[nargs-3]
coords = sys.argv[nargs-2]
print("coordinates file: %s"%coords)

#create links to potentials files
dstO = 'pseudo.O_ONCV_PBE_SG15'
srcO = sys.argv[-1] + '/' + dstO

dstH = 'pseudo.H_ONCV_PBE_SG15'
srcH = sys.argv[-1] + '/' + dstH

if not os.path.exists(dstO):
  print("Create link to %s"%dstO)
  os.symlink(srcO, dstO)
if not os.path.exists(dstH):
  print("Create link to %s"%dstH)
  os.symlink(srcH, dstH)

#run quench
command = "{} {} -c {} -i {}".format(mpicmd,exe,inp1,coords)
print("Run command: {}".format(command))
output1 = subprocess.check_output(command,shell=True)
lines=output1.split(b'\n')

#analyse output of quench
for line in lines:
  if line.count(b'##') and line.count(b'O1'):
    print(line)
    words=line.split()
    x0 = eval(words[3])
    y0 = eval(words[4])
    z0 = eval(words[5])

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

print("Check atom is locked...")
tol = 1.e-4
energy = 0.
count = 0
for line in lines:
  if line.count(b'##') and line.count(b'O1'):
    count = count + 1
    print(line)
    words=line.split()
    x = eval(words[3])
    y = eval(words[4])
    z = eval(words[5])
    print("O1 coordinates = {},{},{}".format(x,y,z))

    if abs(x-x0) > tol or abs(y-y0) > tol or abs(z-z0) > tol:
      print("(x0,y0,z0) = {},{},{}".format(x0,y0,z0))
      print("(x,y,z) = {},{},{}".format(x,y,z))
      sys.exit(1)

if count<3:
  print("O1 coordinates not found in MD output")
  sys.exit(1)

sys.exit(0)
