#!/usr/bin/env python
import sys
import os
import subprocess
import string

print("Test Spin O2 LDA...")

nargs=len(sys.argv)

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
for i in range(4,nargs-5):
  mpicmd = mpicmd + " "+sys.argv[i]
print("MPI run command: {}".format(mpicmd))

exe = sys.argv[nargs-4]
inp = sys.argv[nargs-3]
coords = sys.argv[nargs-2]
print("coordinates file: %s"%coords)

#create links to potentials files
dst1 = 'pseudo.O_ONCVPSP_LDA'
src1 = sys.argv[-1] + '/' + dst1

if not os.path.exists(dst1):
  print("Create link to %s"%dst1)
  os.symlink(src1, dst1)

#run mgmol
command = "{} {} -c {} -i {}".format(mpicmd,exe,inp,coords)
print("Run command: {}".format(command))
output = subprocess.check_output(command,shell=True)

#analyse mgmol standard output

lines=output.split(b'\n')

#check converged energy
energy = 0.
for line in lines:
  if line.count(b'%%'):
    print(line)
    words=line.split()
    energy = eval(words[5][:-1])

ref_energy = -31.6130
print("energy = {}".format(energy))
if abs(ref_energy-energy) > 1.e-3:
  print("Incorrect energy, expected {}".format(ref_energy))
  sys.exit(1)

#make sure forces are below tolerance
tol = 4.e-4
Fz  = 1.06e-2
for line in lines:
  #find output lines with forces
  if line.count(b'##'):
    words=line.split()
    if len(words)==8:
      print(line)
      #check forces in x, y directions below tolerance
      for i in range(5,7):
        force = eval(words[i])
        if abs(force)>tol:
          print("force = {}".format(force))
          sys.exit(1)
      #check value of force in z direction
      if abs(eval(words[7])-Fz)>tol:
          print("force in z dir = {}".format(eval(words[7])))
          sys.exit(1)
      else:
        #switch sign for next test
        Fz = -Fz

