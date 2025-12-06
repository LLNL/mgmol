#!/usr/bin/env python
import sys
import os
import subprocess
import string

print("Test Cl2_ONCVPSP_LDA...")

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
dst1 = 'pseudo.Cl_ONCVPSP_LDA'
src1 = sys.argv[-1] + '/' + dst1

if not os.path.exists(dst1):
  print("Create link to %s"%dst1)
  os.symlink(src1, dst1)

#run mgmol
command = "{} {} -c {} -i {}".format(mpicmd,exe,inp,coords)
print("Run command: {}".format(command))
output = subprocess.check_output(command,shell=True)

#analyse mgmol standard output
#make sure forces are below tolerance
lines=output.split(b'\n')

tol = 4.e-6
Fz  = -7.33e-04
for line in lines:
  num_matches = line.count(b'%%')
  if num_matches:
    print(line)
  #find output lines with forces
  num_matches = line.count(b'##')
  if num_matches:
    words=line.split()
    if len(words)==8:
      print(line)
      #check forces in x, y directions below tolerance
      for i in range(5,7):
        force = eval(words[i])
        if abs(force)>tol:
          print("Force larger than tol, force = {}".format(force))
          sys.exit(1)
      #check value of force in z direction
      if abs(eval(words[7])-Fz)>2.e-5:
          print("force in z dir = {}".format(eval(words[7])))
          sys.exit(1)
      else:
        #switch sign for next test
        Fz = -Fz

