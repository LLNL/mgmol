#!/usr/bin/env python
import sys
import os
import subprocess
import string

print("Test SiH4...")

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
dst1 = 'pseudo.Si'
dst2 = 'pseudo.H'
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
#make sure forces are below tolerance
lines=output.split(b'\n')

tol = 5.e-3
found_forces = False
for line in lines:
  if line.count(b'%%'):
    print(line)
  if line.count(b'##'):
    words=line.split()
    if len(words)==8:
      print(line)
      found_forces = True
      for i in range(5,8):
        if abs(eval(words[i]))>tol:
          sys.exit(1)

if (not found_forces):
  print("no forces found")
  sys.exit(1)

sys.exit(0)
