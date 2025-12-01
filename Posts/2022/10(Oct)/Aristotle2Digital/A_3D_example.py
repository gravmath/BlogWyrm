import matplotlib.pyplot as plt
import numpy             as np
import sympy             as sp


#example from https://sites.calvin.edu/scofield/courses/m355/handouts/definiteMatrices.pdf
A      = np.array([[1,1,0],[1,2.5,1],[0,1,3.5]])
print(A - A.transpose())
print(np.linalg.eigvals(A))
print(np.linalg.det(A))
print(np.pi**(1.5)/np.sqrt(np.linalg.det(A)))

x, y, z = sp.symbols('x y z')
r          = sp.Matrix([[x],[y],[z]])
quad_form  = sp.simplify((sp.transpose(r)*A*r)[0,0])

print(quad_form)

def my_exp(x,y,z):
    arg = 1.0*x**2 + 2.0*x*y + 2.5*y**2 + 2.0*y*z + 3.5*z**2
    return np.exp(-arg)

N      = 100_000_00
mins   = {'x':-6.0,'y':-24.0,'z':-24.0,}
maxs   = {'x':6.0, 'y':24.0, 'z':24.0,}
deltas = {}
vol    = 1
for k in mins.keys():
    deltas[k] = maxs[k] - mins[k]
    vol      *= deltas[k]
#import pdb; pdb.set_trace()

X = mins['x'] + deltas['x']*np.random.rand(N)
Y = mins['y'] + deltas['y']*np.random.rand(N)
Z = mins['z'] + deltas['z']*np.random.rand(N)

I0 = my_exp(X,Y,Z)*vol

print(np.mean(I0))
