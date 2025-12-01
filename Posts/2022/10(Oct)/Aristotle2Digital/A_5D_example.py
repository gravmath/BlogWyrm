import matplotlib.pyplot as plt
import numpy             as np


def my_exp(v,w,x,y,z):
    return np.exp(-(v*v + x*y + 2*y*y))

A      = np.array([[20,0,0,1,1/2],[0,2,0,0,2],[0,0,1,2,0],[1,0,2,4,1],[1/2,2,0,1,4]])
print(A - A.transpose())
print(np.linalg.eigvals(A))

'''
N      = 10_000
mins   = {'v':,'w':,'x':-6.0,'y':-12.0,'z':,}
maxs   = {'x':6.0, 'y':12.0}
deltas = {}
vol    = 1
for k in mins.keys():
    deltas[k] = maxs[k] - mins[k]
    vol      += deltas[k]


X = mins['x'] + deltas['x']*np.random.rand(N)
Y = mins['y'] + deltas['y']*np.random.rand(N)

I0 = my_exp(X,Y)*vol
I  = X**4*Y**2*I0*vol

print(np.mean(I0),2.0*np.pi/np.sqrt(7.0))
print(np.mean(I/I0),144.0/343.0)

X_anal = np.arange(-2,2.1,0.1)
Y_anal = np.arange(-4,4.1,0.1)

Z_anal = np.zeros((len(X_anal),len(Y_anal)))
for i,x in enumerate(X_anal):
    for j,y in enumerate(Y_anal):
        Z_anal[i,j] = (x**1)*my_exp(x,y)
print(np.max(Z_anal))

#import pdb; pdb.set_trace()
plt.contour(X_anal,Y_anal,Z_anal.T,levels=10)
plt.xlabel('x')
plt.ylabel('y')
plt.show()        

N  = np.array([[1,1/2],[1/2,2]])
Vs = np.linalg.eigvals(N)

print(Vs)
'''