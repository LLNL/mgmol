# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
# Python program to generate input coordinates for mgmol from mgmol output
#
# use: python mgmol2lrs.py mgmol_output > lrs.in
#-------------------------------------------------------------------------------
import sys, string
import outputTools

output=open(sys.argv[1],'r')
Lo    =output.readlines()
lo=len(Lo)  ## no of lines in file

searchterm1='centers'
searchterm2='spreads'

for line in range(lo): ## loop over lines of file 
  num_matches1 = string.count(Lo[line], searchterm1)
  num_matches2 = string.count(Lo[line], searchterm2)
  if num_matches2 and num_matches1:
      lrs=[]
      flag1=0
      flag2=0
      #print 'loop starting at', line+1
      for line2 in range(line+1,lo):
        words=string.split(Lo[line2])
        if len(words)>0:
          if words[0]=='&&': 
            flag1=1
            flag2=1
            x=words[-4]
            y=words[-3]
            z=words[-2]
            #print name+'\t'+x+'\t'+y+'\t'+z
            lrs.append(x.ljust(10)+y.ljust(10)+z.ljust(10))
          else:
            flag2=0
          if flag1!=flag2: break

#print new lrs file
for lr in lrs:
  print lr

