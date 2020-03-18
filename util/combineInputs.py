# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
import sys, string
from re import split

input1=open(sys.argv[1],'r')
input2=open(sys.argv[2],'r')
axis=eval(sys.argv[3])

n=[1,1,1]
n[axis]=2

L1=input1.readlines()
L2=input2.readlines()

# list possible species in a dictionary with their count
species={'H':0,'He':0,'Li':0,'Be':0,'B':0,'C':0,'N':0,'O':0,'F':0,'Na':0,
         'Mg':0,'Al':0,'Si':0,'P':0,'S':0,'Cl':0,'K':0,'Ca':0,'D':0,'In':0}

numbers_list=['0','1','2','3','4','5','6','7','8','9']

########################################
def get_sp(name):
  name1=name[0:1]
  name2=name[0:2]
  sp = 'none'
  if name1 in species.keys():
    sp = name1
  if name2 in species.keys():
    sp = name2

  return sp

########################################
def next_line(L,line,print_flag,check_ret):
  if line<0:
    return -1
  next_line=line
  max_nlines=len(L)
  while True:
    if next_line>=max_nlines:
      return -1
    new_line=L[next_line]
    if len(new_line)>0:
      first_char=new_line[0:1]
      if first_char!='#':
        if not check_ret:
          break
        elif first_char!='\n' and first_char!=' ':
          #print 'First character:',first_char,' in line of length ',len(new_line)
          break
    if print_flag:
      print L[next_line],
    next_line=next_line+1
    
  return next_line
  
########################################

pseudopotential_list=[]
pseudopotential_list1=[]
pseudopotential_list2=[]

for line in L1:
  words=string.split(line)
  if len(words)>0:
    word=string.split(words[0],".")
    if word[0]=='pseudo':
      #print words[0],word,'pseudo'
      pseudopotential_list1.append(words[0])
      pseudopotential_list.append(words[0])

for line in L2:
  words=string.split(line)
  if len(words)>0:
    word=string.split(words[0],".")
    if word[0]=='pseudo':
      #print words[0],word,'pseudo'
      pseudopotential_list2.append(words[0])
      if words[0] not in pseudopotential_list:
        pseudopotential_list.append(words[0])

#print pseudopotential_list
#print pseudopotential_list1
#print pseudopotential_list2

######################################################################
#establish species number for new (global) list of species
species_map1={}
species_map2={}

index=1
for pot in pseudopotential_list1:
  species_map1[index]=pseudopotential_list.index(pot)+1
  index=index+1
  
index=1
for pot in pseudopotential_list2:
  species_map2[index]=pseudopotential_list.index(pot)+1
  index=index+1

#print species_map1
#print species_map2

######################################################################
count=0
count_atom = 0
count_lrs = 0
reached_atoms = 0
nsp = 0
l=len(L1)  ## no of lines in file



  
line1=0
line2=0
nst=-1
check_for_ret=True
for line in L1: ## loop over lines of file
  #print '####################'
  line1=next_line(L1,line1,True,check_for_ret)
  line2=next_line(L2,line2,False,check_for_ret)
  #print '1:',L1[line1]
  #print '2:',L2[line2]
  
  if line1<0:
    word1=[]
  else:
    word1=string.split(L1[line1])
  
  if line2<0:
    word2=[]
  else:
    word2=string.split(L2[line2])
  
  if line1<0 and line2<0:
    break
  
  count=count+1
  #print 'count=',count
  #print 'word1=',word1
  #print 'word2=',word2

  if count==3:
    nnx=n[0]*eval(word1[0])
    nny=n[1]*eval(word1[1])
    nnz=n[2]*eval(word1[2])
    print nnx,'\t',nny,'\t',nnz
  elif  count==2:
    b1=eval(word1[0])
    b2=eval(word1[1])
    b3=eval(word1[2])
    e1=eval(word1[3])
    e2=eval(word1[4])
    e3=eval(word1[5])
    d1=e1-b1
    d2=e2-b2
    d3=e3-b3
    print b1,'\t',b2,'\t',b3,'\t',
    print e1+(n[0]-1)*d1,'\t',e2+(n[1]-1)*d2,'\t',e3+(n[2]-1)*d3
  elif  count==4:
    print len(pseudopotential_list)
  elif  count==5:
    for pot in pseudopotential_list:
      print pot,' 1'
    print ' '
    line1=line1+len(pseudopotential_list1)
    line2=line2+len(pseudopotential_list2)
  elif  count==33:
    ne1=eval(word1[0])
    ne2=eval(word2[0])
    nst=ne1+ne2
    print nst
  elif  count==34:
    ne1=eval(word1[0])
    ne2=eval(word2[0])
    print ne1+ne2,'  ',0.5*(eval(word1[1])+eval(word2[1]))
  elif  count==35:
    ni1=eval(word1[0])
    ni2=eval(word2[0])
    print ni1+ni2
  elif count>41 and len(word1)>4:
    check_for_ret=False
    myspecies1 = get_sp(word1[0])
    myspecies2 = get_sp(word2[0])
    if myspecies1 !='none':
      reached_atoms = 1
      
      x=eval(word1[2])
      y=eval(word1[3])
      z=eval(word1[4])
      count_atom = count_atom + 1
      name = myspecies1 + str(count_atom)
      sp=species_map1[ eval(word1[1]) ]
      print name.ljust(7),str(sp).rjust(3),str(x).rjust(16),str(y).rjust(16),str(z).rjust(16),
      if len(word1)>5:
        for word in word1[5:]:
          print word.rjust(3),
      print  '\n',
    else:
      print L1[line1],  
      
    if myspecies2 !='none':
      reached_atoms = 1
      
      x=eval(word2[2])+(n[0]-1)*d1
      y=eval(word2[3])+(n[1]-1)*d2
      z=eval(word2[4])+(n[2]-1)*d3
      count_atom = count_atom + 1
      name = myspecies2 + str(count_atom)
      sp=species_map2[ eval(word2[1]) ]
      print name.ljust(7),str(sp).rjust(3),str(x).rjust(16),str(y).rjust(16),str(z).rjust(16),
      if len(word2)>5:
        for word in word2[5:]:
          print word.rjust(3),
      print  '\n',
    else:
      print L2[line],  

  else:
    #localization centers
    if ( len(word1)==4 or len(word1)==3 ) and reached_atoms > 0:
      count_lrs=count_lrs+1
      x=eval(word1[0])
      y=eval(word1[1])
      z=eval(word1[2])
      if len(word1)==4:
        print str(x).ljust(16),str(y).ljust(16),str(z).ljust(16),word[3].rjust(4)
      else:
        print str(x).ljust(16),str(y).ljust(16),str(z).ljust(16)
    elif len(word1)>0:
      print L1[line1],
    
    if ( len(word2)==4 or len(word2)==3 ) and reached_atoms > 0:
      count_lrs=count_lrs+1
      #print word1,word2,count_lrs
      x=eval(word2[0])+(n[0]-1)*d1
      y=eval(word2[1])+(n[1]-1)*d2
      z=eval(word2[2])+(n[2]-1)*d3
      if len(word2)==4:
        print str(x).ljust(16),str(y).ljust(16),str(z).ljust(16),word[3].rjust(4)
      else:
        print str(x).ljust(16),str(y).ljust(16),str(z).ljust(16)
    #print 'lrs: ',count_lrs
  line1=line1+1
  line2=line2+1
  
  if count_lrs==nst:
    break
