#!/usr/bin/env python
import sys
import os
import subprocess
import string

print("Test SiH4...")

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]
coords = sys.argv[6]
print("coordinates file: %s"%coords)

#create links to potentials files
dst1 = 'pseudo.Si'
dst2 = 'pseudo.H'
src1 = sys.argv[7] + '/' + dst1
src2 = sys.argv[7] + '/' + dst2

if not os.path.exists(dst1):
  print("Create link to %s"%dst1)
  os.symlink(src1, dst1)
if not os.path.exists(dst2):
  print("Create link to %s"%dst2)
  os.symlink(src2, dst2)

#run mgmol
command = "{} {} -c {} -i {}".format(mpicmd,exe,inp,coords)
output = subprocess.check_output(command,shell=True)

#analyse mgmol standard output
#make sure forces are below tolerance
lines=output.split(b'\n')

tol = 5.e-3
for line in lines:
  num_matches = line.count(b'%%')
  if num_matches:
    print(line)
  num_matches = line.count(b'##')
  if num_matches:
    words=line.split()
    if len(words)==8:
      print(line)
      for i in range(5,8):
        if abs(eval(words[i]))>tol:
          sys.exit(1)

