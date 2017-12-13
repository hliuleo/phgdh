from researchcode.plotting.plot_set import *
import sys
from numpy.linalg import norm
from numpy import dot
import math

def angle_between(a,b):
  arccosInput = dot(a,b)/norm(a)/norm(b)
  arccosInput = 1.0 if arccosInput > 1.0 else arccosInput
  arccosInput = -1.0 if arccosInput < -1.0 else arccosInput
  return math.acos(arccosInput)/math.pi*180

coms = np.loadtxt(sys.argv[1])
coms = coms.reshape((len(coms), 5, 3))

def calAngles(d):
    def getNV(p):
        v1 = p[1]-p[0]
        v2 = p[2]-p[0]
        nv = np.cross(v1, v2)
        return nv
    tp1 = d[[0,2,3]]
    tp2 = d[[1,2,3]]
    op1 = d[[0,2,4]]
    op2 = d[[1,2,4]]
    tangles = angle_between(getNV(tp1), getNV(tp2))
    oangles = angle_between(getNV(op1), getNV(op2))
    return [oangles, tangles]

coords = []
for c in coms:
    coords.append(calAngles(c))
np.savetxt('coords.dat', coords)
