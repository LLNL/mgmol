# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
#
#usage: python generate_LinXYZPot.py -b -f cell.dat -e 0.001 -d 0 -o outputfile.xyz
#
import sys, string, struct, getopt
from array import array

opts, args = getopt.getopt(sys.argv[1:], "hf:e:d:bo:")

def usage():
  usage = """
  -h --help                  Prints this
  -f --file (argument)       input file
  -e --efield (argument)     field
  -d --direction (argument)  Direction
  -o --output (argument)     Output
  -b --binary                binary output
  """
  print usage

bin_flag=0
ifilename='0'
ofilename='0'
for o,a in opts:
  #print o,a
  if o in ("-h", "--help"):
    usage()
    sys.exit()
  elif o in ("-f", "--file"):
    ifilename = a
    print "Input file is", ifilename
  elif o in ("-e", "--efield"):
    efield = eval(a)
    print "Field is" , efield
  elif o in ("-d", "--direction"):
    direction = eval(a)
    print "Direction is" , direction
  elif o in ("-o", "--output"):
    print "Output file is " , a
    ofilename=a
  elif o in ("-b", "--binary"):
    bin_flag=1
  else:
    print 'unknown option'
    print opts
    usage()
    sys.exit()

ifile      =open(ifilename,'r')

if bin_flag>0:
  print 'Binary output...'
  output=open(ofilename,'wb')
else:
  output=open(ofilename,'w')

flag=0
origin=[0.,0.,0.]
end=[0.,0.,0.]
while flag<1:
  line = ifile.readline()
  words=string.split(line)
  if words[0][0]!='#':
    origin[0]=eval(words[0])
    origin[1]=eval(words[1])
    origin[2]=eval(words[2])
    end[0]=eval(words[3])
    end[1]=eval(words[4])
    end[2]=eval(words[5])
    print 'origin=',origin
    print 'end=',end
    flag=1

if bin_flag>0:
  packed_string = struct.pack('fff',*origin)
  output.write(packed_string)
  packed_string = struct.pack('fff',*end)
  output.write(packed_string)
  #float_array = array('f', [0.,-4.69,-14.07])
  #float_array = array('f', origin)
  #float_array.tofile(output)
else:
  output.write('#Cell\n')
  output.write(str(origin[0]))
  output.write('\t')
  output.write(str(origin[1]))
  output.write('\t')
  output.write(str(origin[2]))
  output.write('\t')
  output.write(str(end[0]))
  output.write('\t')
  output.write(str(end[1]))
  output.write('\t')
  output.write(str(end[2]))
  output.write('\n')

flag=0
n=[0,0,0]
while flag<1:
  line = ifile.readline()
  words=string.split(line)
  if words[0][0]!='#':
    n[0]=eval(words[0])
    n[1]=eval(words[1])
    n[2]=eval(words[2])
    print 'mesh=',n
    flag=1

if bin_flag>0:
  packed_string = struct.pack('iii',*n)
  output.write(packed_string)
else:
  output.write('#Mesh:\n')
  output.write(str(n[0]))
  output.write('\t')
  output.write(str(n[1]))
  output.write('\t')
  output.write(str(n[2]))
  output.write('\n')
  output.write('#E-field=')
  output.write(str(efield))
  output.write(', direction=')
  output.write(str(direction))
  output.write('\n')

h=((end[0]-origin[0])/n[0],(end[1]-origin[1])/n[1],(end[2]-origin[2])/n[2])

count=0
for i in range(n[0]):
  for j in range(n[1]):
    for k in range(n[2]):
      
      val=0.
      if direction==0:
        val=(origin[direction]+i*h[direction])*efield
      if direction==1:
        val=(origin[direction]+j*h[direction])*efield
      if direction==2:
        val=(origin[direction]+k*h[direction])*efield
      if bin_flag>0:
        packed_string = struct.pack('f',val)
        output.write(packed_string)
      else:
        output.write(str(val))
        output.write('\n')
      count=count+1

if count!=(n[0]*n[1]*n[2]):
  print 'ERROR: count!=n[0]*n[1]*n[2]'
