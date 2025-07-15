#!/usr/bin/env python
import sys
import os
import subprocess
import string
import shutil

print("Test Override coordinates...")

nargs=len(sys.argv)

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
for i in range(4,nargs-6):
  mpicmd = mpicmd + " "+sys.argv[i]
print("MPI run command: {}".format(mpicmd))

exe = sys.argv[nargs-6]
inp1 = sys.argv[nargs-5]
inp2 = sys.argv[nargs-4]
coords1 = sys.argv[nargs-3]
coords2 = sys.argv[nargs-2]
print("coordinates file: %s"%coords1)

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

#run quench
command = "{} {} -c {} -i {}".format(mpicmd,exe,inp1,coords1)
print("Run command: {}".format(command))
output1 = subprocess.check_output(command,shell=True)
lines=output1.split(b'\n')

#analyse output of quench
for line in lines:
  if line.count(b'%%'):
    print(line)

#run quench with shifted coordinates
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

command = "{} {} -c {} -i {}".format(mpicmd,exe,inp2,coords2)
print(command)
output2 = subprocess.check_output(command,shell=True)

#remove created files
shutil.rmtree(restart_file)
os.remove('wave.out')

#analyse mgmol standard output
lines=output2.split(b'\n')

flag = 0
for line in lines:
  if line.count(b'%%'):
    print(line)
  if line.count(b'achieved'):
    flag=1

if flag==0:
  print("second run did not converge...")
  sys.exit(1)

sys.exit(0)
