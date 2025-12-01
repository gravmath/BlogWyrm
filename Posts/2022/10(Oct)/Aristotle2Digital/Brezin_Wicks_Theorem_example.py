import matplotlib.pyplot as plt
import numpy             as np


def my_exp(x,y):
    return np.exp(-(x*x + x*y + 2*y*y))

N      = 10_000_000
mins   = {'x':-6.0,'y':-6.0}
maxs   = {'x':6.0, 'y':6.0}
deltas = {'x':maxs['x']-mins['x'],
          'y':maxs['y']-mins['y'],}
vol    = deltas['x']*deltas['y']

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