# Generate local potential according to
# Goedecker, Teter, Hutter, Phys. rev. B 54 (3), 1996
# Parameters from M. Krack, Theor. Chem. Acc. 114, 2005
from math import exp, erf, sqrt, pi

#coefficients for H
rloc = 0.4
c1 = -14.081155
c2 =   9.626220
c3 =  -1.783616
c4 =   0.085152
zion = 3.
anumber = 3
name = "LithiumGTH_PBE"
mass = 3.

def radialfunction(r):
  alpha = (r/rloc)**2
  val = exp(-0.5*alpha)*(c1+c2*alpha+c3*alpha*alpha+c4*alpha*alpha*alpha)
  if r>1.e-8:
    val = val - zion*erf(r/(sqrt(2.)*rloc))/r
  else:
    #print("special case for r = {}".format(r))
    val = val -zion*sqrt(2.)/(sqrt(pi)*rloc)

  return val

npts = 301

#header
print("# Short description of the species type. One line only!")
print(name)
print("#")
print("White")
print("#radii of balls and covalent bonds")
print("0.4 1.0")
print("# Nlcc flag")
print("0")
print("# Atomic number")
print(anumber)
print("# Atomic mass")
print(mass)
print("# Number of valence electrons")
print(zion)
print("# Gaussian core charge parameter rc")
print("1.")
print("# Number of potentials")
print("1")
print("# l-value for state which is local")
print("0 0")
print("# Local potential radius")
print("3.")
print("# Non-local potential radius")
print("3.")
print("# number of points in radial grid")
print(npts)
print("# log mesh parameter")
print("0.")
print("# radial grid, reference state, and potential for l=0")

#potential
for i in range(npts):
  r = round(0.01*i,2)
  f = radialfunction(r)
  print("{} {}".format(r,f))
