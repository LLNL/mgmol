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

for line in lines1:
  if line.count('Timer:'):
    words=line.split()
    timers1[words[1]]=words[6]

for line in lines2:
  if line.count('Timer:'):
    words=line.split()
    timers2[words[1]]=words[6]

#analyse timers
results={}
for timer in timers1.keys():
  if timer in timers2.keys():
    #compute relative difference
    diff=(eval(timers2[timer])-eval(timers1[timer]))/eval(timers1[timer])
    if abs(diff)>rel_threshold and abs(eval(timers1[timer]))>abs_threshold:
      results[timer]=diff

sorted_timers=sorted(results.items(), key=operator.itemgetter(1))

#print results
print('-----------------------------------------------------------------------------')
print('Timer                                     time1    time2    relative diff.(%)')
print('-----------------------------------------------------------------------------')
for timer in reversed(sorted_timers):
  key=timer[0]
  print(key.ljust(50)),
  print(str(eval(timers1[key])).ljust(10)),
  print(str(eval(timers2[key])).ljust(10)),
  print(str(100.*timer[1]).ljust(20))
