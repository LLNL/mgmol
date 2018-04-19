# use: python mgmol2input.py coords.in  dz > coords_new.in
#-------------------------------------------------------------------------------
import sys, string

input1 = open(sys.argv[1],'r')
dz     = eval(sys.argv[2])

#print new input file
Li=input1.readlines()

for line in Li: ## loop over lines of file
  words=string.split(line)
  if len(words)>0:
    name=words[0]
    if name=="Cl2" and len(words)>4:
      z=eval(words[4])+dz
      print words[0]+'\t'+words[1]+'\t'+words[2]+'\t'+words[3]+'\t'+str(z)
    else:
      print line,
  else:
    print line,
