import matplotlib.pyplot as plt
import numpy             as np
import sympy             as sp


#example from https://sites.calvin.edu/scofield/courses/m355/handouts/definiteMatrices.pdf
A      = np.array([[3,1,5,3],[1,5,2,0],[5,2,10,3],[3,0,3,14]])
#print(A - A.transpose())
print(np.linalg.eigvals(A))
print(np.pi**2/np.sqrt(np.linalg.det(A)))

w, x, y, z = sp.symbols('w x y z')
r          = sp.Matrix([[w],[x],[y],[z]])
quad_form  = sp.simplify((sp.transpose(r)*A*r)[0,0])

#print(quad_form)

def my_exp(w,x,y,z):
    arg = 3*w**2 + 2*w*x + 10*w*y + 6*w*z + 5*x**2 + 4*x*y + 10*y**2 + 6*y*z + 14*z**2
    return np.exp(-arg)


N      = 100_000_000
mins   = {'w':-20.0,'x':-20.0,'y':-20.0,'z':-20.0,}
maxs   = {'w':20.0, 'x':20.0, 'y':20.0, 'z':20.0,}
deltas = {}
vol    = 1
for k in mins.keys():
    deltas[k] = maxs[k] - mins[k]
    vol      *= deltas[k]
#import pdb; pdb.set_trace()

W = mins['w'] + deltas['w']*np.random.rand(N)
X = mins['x'] + deltas['x']*np.random.rand(N)
Y = mins['y'] + deltas['y']*np.random.rand(N)
Z = mins['z'] + deltas['z']*np.random.rand(N)

I0 = my_exp(W,X,Y,Z)*vol

print(np.mean(I0))


