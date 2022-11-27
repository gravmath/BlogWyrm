import matplotlib.pyplot as plt
import numpy             as np
import sympy as sym

bx, by, bz, Bxx, Bxy, Bxz, Byy, Byz, Bzz = sym.symbols('bx by bz Bxx Bxy Bxz Byy Byz Bzz')

expon = ( bx*bx*Bxx + 2*bx*by*Bxy + 2*bx*bz*Bxz + by*by*Byy + 2*by*bz*Byz + bz*bz*Bzz )/2
dist  = sym.exp(expon)

Ex2yz3 = sym.diff(sym.diff(sym.diff(dist,bx,2),by,1),bz,3)
print('E[x^2yz^3]: ',Ex2yz3.subs({bx:0,by:0,bz:0}))