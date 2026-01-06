# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#usage:
#python compareTimers.py output1 output2

import sys, string, operator

#ignore timers smaller than abs_threshold
abs_threshold=5.

#ignore timers of relative difference less than rel_threshold
rel_threshold=0.05

#read files
input1=open(sys.argv[1],'r')
input2=open(sys.argv[2],'r')

lines1=input1.readlines()
lines2=input2.readlines()

#extract timers from files
timers1={}
timers2={}

count1=0
count2=0

for line in lines1:
  if line.count('comp_res_from_Hphi'):
    words=line.split()
    count1=eval(words[8])

print("count 1 = {}".format(count1))

for line in lines2:
  if line.count('comp_res_from_Hphi'):
    words=line.split()
    count2=eval(words[8])

print("count 2 = {}".format(count2))

for line in lines1:
  if line.count('Timer:'):
    words=line.split()
    key=words[1]
    val=eval(words[6])
    if 'gemm' in key:
      key='gemm'
    timers1[key]=val/count1

for line in lines2:
  if line.count('Timer:'):
    words=line.split()
    key=words[1]
    val=eval(words[6])
    if 'gemm' in key:
      key='gemm'
    timers2[key]=val/count2

#analyse timers
results={}
for key in timers1.keys():
  if key in timers2.keys():
    #compute relative difference
    diff=(timers2[key]-timers1[key])/timers1[key]
    if abs(diff)>rel_threshold and abs(timers1[key])>abs_threshold/count1:
      results[key]=diff
    else:
      if 'MGmol::total' in key:
        results[key]=diff
      if 'quench' in key:
        results[key]=diff
      if 'gemm' in key:
        results[key]=diff

sorted_timers=sorted(results.items(), key=operator.itemgetter(1))

#print results
print('---------------------------------------------------------------------------------------')
print('Timer                                             time1     time2     relative diff.(%)')
print('---------------------------------------------------------------------------------------')
ndec=2
for timer in reversed(sorted_timers):
  key=timer[0]
  print(key.ljust(50), end="")
  print(str(round(timers1[key]*count1,ndec))[:8].ljust(10), end="")
  print(str(round(timers2[key]*count1,ndec))[:8].ljust(10), end="")
  print(str(round(100.*timer[1],ndec))[:8].ljust(10))
