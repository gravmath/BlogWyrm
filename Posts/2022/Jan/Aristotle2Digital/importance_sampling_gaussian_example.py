import matplotlib      as mpl
import numpy           as np
import pandas          as pd 
import scipy.special   as special
import scipy.integrate as integ

import matplotlib.pyplot as plt

Ns    = [i for i in range(10,2010,10)]
xmin  = 0
xmax  = 10
df    = pd.DataFrame()
exact = np.sqrt(np.pi)/2.0*special.erf(10)  #exact answeer for this integral in terms of the error function

def f(x):     return np.exp(-x*x)
def pu(span): return 1/(span[1]-span[0])
def pe(x):    return np.exp(-x)
def xk(x):    return -np.log(1-x)

for N in Ns:
    U  = np.random.random((N,))
    FU = f(xmin + xmax*U)
    PU = pu((xmin,xmax))
    XK = xk(U)
    PE = pe(XK)
    FX = f(XK)
    
    mean_U = np.mean(FU/PU)
    sig_U  = np.std(FU/PU)/np.sqrt(N)
    err_U  = mean_U - exact
    mean_K = np.mean(FX/PE)
    sig_K  = np.std(FX/PE)/np.sqrt(N)
    err_K  = mean_K - exact

    row = {'N':N,'I_U':mean_U,'sig_I_U':sig_U,'error_U':err_U,'I_K':mean_K,'sig_I_K':sig_K,'error_K':err_K}
    df  = df.append(row,ignore_index=True)

X = np.linspace(xmin,xmax,100)
F = f(X)

print('The "true" answer is:                  ',np.sqrt(np.pi)/2.0*special.erf(10))
print('The numerical integration answer is:   ',integ.quad(f,0,10)[0])
print(df)

plt.plot(X,F,'r-')
plt.plot(X,np.exp(-X),'k-')
plt.xlabel('$x_{y^2}$')
plt.show()

fig_h = plt.figure()
ax_h  = fig_h.add_subplot(111)
ax_h.hist(XK,color='#ff6666')
ax_h.set_xlabel('$X$',fontsize=16)
ax_h.set_ylabel('Number of Hits $N$',fontsize=16)
#fig_h.savefig('C:/Users/byecs/Documents/GitHub/BlogWyrm/Posts/2022/Jan/Aristotle2Digital/exponential_dist.png')
print('mean of x:',np.mean(XK),' with a standard deviation of: ',np.std(XK),XK.shape)

fig_err_by_N = plt.figure()
ax_err_by_N  = fig_err_by_N.add_subplot(111)
ax_err_by_N.plot(Ns,df['error_U'],'ko-',label='uniform')
ax_err_by_N.plot(Ns,df['error_K'],'rs-',label='exponential')
ax_err_by_N.set_xlabel('$N$',fontsize=16)
ax_err_by_N.set_ylabel('$Error \;\; (\epsilon$)',fontsize=16)
ax_err_by_N.legend()
#fig_err_by_N.savefig('C:/Users/byecs/Documents/GitHub/BlogWyrm/Posts/2022/Jan/Aristotle2Digital/error_uniform_v_exponential_sampling.png')
plt.show()

fig_var_by_N = plt.figure()
ax_var_by_N  = fig_var_by_N.add_subplot(111)
ax_var_by_N.plot(Ns,df['sig_I_U'],'ko-',label='uniform')
ax_var_by_N.plot(Ns,df['sig_I_K'],'rs-',label='exponential')
ax_var_by_N.set_xlabel('$N$',fontsize=16)
ax_var_by_N.set_ylabel('$Standard \;\; Deviation \;\; (\sigma$)',fontsize=16)
ax_var_by_N.legend()
#fig_var_by_N.savefig('C:/Users/byecs/Documents/GitHub/BlogWyrm/Posts/2022/Jan/Aristotle2Digital/var_uniform_v_exponential_sampling.png')
plt.show()