#!/usr/bin/env python
import sys
import os
import subprocess
import string

def check_is_small(f):
  tol = 3.e-5
  print(f)
  if abs(eval(f))>tol:
    print("value not small")
    return 1
  return 0

def check_opposed(f0,f1):
  tol = 1.e-5
  print("check opposites: {} {}".format(f0,f1))
  if abs(eval(f0)+eval(f1))>tol:
    print("valuee not opposites")
    return 1
  return 0

def check_equal(f0,f1):
  tol = 1.e-5
  print("check equal : {} {}".format(f0,f1))
  if abs(eval(f0)-eval(f1))>tol:
    print("valuee not equal")
    return 1
  return 0

#################################################################
print("Test Forces...")

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
dst = 'pseudo.P_ONCV_PBE_SG15'
src = sys.argv[nargs-1] + '/' + dst

if not os.path.exists(dst):
  print("Create link to %s"%dst)
  os.symlink(src, dst)

#run mgmol
command = "{} {} -c {} -i {}".format(mpicmd,exe,inp,coords)
print("Run command: {}".format(command))
output = subprocess.check_output(command,shell=True)

os.remove(dst)

#analyse mgmol standard output
lines=output.split(b'\n')

forces=[]
found_forces = False
for line in lines:
  if line.count(b'%%'):
    print(line)
  if line.count(b'##'):
    words=line.split()
    if len(words)==8:
      print(line)
      found_forces = True
      f=[words[5],words[6],words[7]]
      forces.append(f)

if (not found_forces):
  print("no forces found")
  sys.exit(1)

if len(forces)<4:
  print("Needs 4 forces")
  sys.exit(1)

print("Check values of forces...")
#some forces should be 0 by symmetry
ret = check_is_small(forces[2][0])
ret = ret + check_is_small(forces[3][0])
ret = ret + check_is_small(forces[0][1])
ret = ret + check_is_small(forces[1][1])

#some forces should be opposed to each other by symmetry
ret = ret + check_opposed(forces[0][0],forces[1][0])
ret = ret + check_opposed(forces[2][1],forces[3][1])
ret = ret + check_opposed(forces[0][2],forces[2][2])

#some forces should be the same by symmetry
ret = ret + check_equal(forces[0][2],forces[1][2])
ret = ret + check_equal(forces[2][2],forces[3][2])

if ret > 0:
  sys.exit(1)

sys.exit(0)
