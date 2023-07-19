# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
import sys, string

def countNumAtoms(ifile):
  ifile.seek(0)
  lines=ifile.readlines()

  # list possible species in a dictionary with their count
  species={'H':0,'D':0,'He':0,'Li':0,'Be':0,'B':0,'C':0,'N':0,'O':0,'F':0,'Na':0,
          'Mg':0,'Al':0,'Si':0,'P':0,'S':0,'Cl':0,'Ca':0,'In':0,'Au':0}
  flag1=0
  flag2=0

  for line in lines: ## loop over lines of file 
    i=0
    for c in line:
      flag1=0
      if c=='#': 
        if line[i+1]=='#':
          if line[i+2]!='#':
            flag1=1
            flag2=1 #flag2 turned on when first '##' found
            ss = line.split()
            one = ss[-7]
            two   = ss[-6]
            three = ss[-5]
            four  = ss[-4]
            ii=0
            if one[0]=='*':
              ii=1
            if one[ii:ii+2] in species.keys():
              species[one[ii:ii+2]]=species[one[ii:ii+2]]+1
            else:
              species[one[ii]]=species[one[ii]]+1
            break   
          else:
            break
      i=i+1
    #stop when first atom was found, but new line does not contain one
    if flag1!=flag2: break

  na=0
  for sp in species.keys(): 
    na=na+species[sp]

  return na

#################################################################
def countNumMLWC(ifile):
  ifile.seek(0)
  lines=ifile.readlines()
  nmlwc=0

  flag1=0
  flag2=0

  for line in lines: ## loop over lines of file 
    i=0
    for c in line:
      flag1=0
      if c=='&': 
        if line[i+1]=='&':
          if line[i+2]!='&':
            nmlwc=nmlwc+1
            flag1=1
            flag2=1
            break   
          else:
            break
      i=i+1
    if flag1!=flag2: break

  return nmlwc
