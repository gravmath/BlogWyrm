import matplotlib.pyplot as plt
import numpy             as np
import sympy as sym

A_lst = [[2,1],[1,4]]

A, bx, by, Bxx, Bxy, Byy = sym.symbols('A bx by Bxx Bxy Byy')

expon = ( bx*bx*Bxx + 2*bx*by*Bxy + by*by*Byy)/2
dist  = sym.exp(expon)

Ex42 = sym.diff(sym.diff(dist,bx,4),by,2)
print('E[x^4y^2]: ',Ex42.subs({bx:0,by:0}))

Ex40 = sym.diff(dist,bx,4)
print('E[x^4]: ',Ex40.subs({bx:0,by:0}))

A = sym.Matrix(A_lst)
print('A: ',A,'B: ',A.inv())

def my_exp(x,y):
    return np.exp(-(x*x + x*y + 2*y*y))

A_arr  = np.array(A_lst)
B      = np.linalg.inv(A_arr)

N      = 10_000_000
mins   = {'x':-6.0,'y':-6.0}
maxs   = {'x':6.0, 'y':6.0}
deltas = {'x':maxs['x']-mins['x'],
          'y':maxs['y']-mins['y'],}
vol    = deltas['x']*deltas['y']

X = mins['x'] + deltas['x']*np.random.rand(N)
Y = mins['y'] + deltas['y']*np.random.rand(N)

I0 = my_exp(X,Y)*vol
IB = X**4*Y**2*my_exp(X,Y)*vol
IP = X**4*my_exp(X,Y)*vol

print(np.mean(I0),2.0*np.pi/np.sqrt(7.0))
print(np.mean(IB)/np.mean(I0),3*B[0,0]**2*B[1,1] + 12*B[0,0]*B[1,0]**2)
print(np.mean(IP)/np.mean(I0),3*B[0,0]**2)