import matplotlib        as mpl
import matplotlib.pyplot as plt
import numpy             as np
import pandas            as pd
import scipy             as sp
import scipy.integrate   as qint

def double_hump(x,x_min,x_max):
    if x <= x_max and x >= x_min:
        return x*np.sin(x) + np.cos(x)
    else:
        return 0.0

def exp(x,x_min,x_max):
    if x >= 0:
        return np.exp(-x)
    else:
        return 0.0

def cos_pow(x,x_min,x_max):
    if x <= x_max and x >= x_min:
        return (np.cos(x))**6/1.3530503614562723
    else:
        return 0.0

def pos_tanh(x,x_min,x_max):
    if x <= x_max and x >= x_min:
        return (np.abs(np.tanh(x)))/4.2770877257340025
    else:
        return 0.0

x_lim    = 2.8
N_trials = 1_000_000
delta    = 7.0
Ds       = (2*np.random.random(N_trials) - 1.0)*delta
Rs       = np.random.random(N_trials)
n_accept = 0
n_burn   = 1_000
x_0      = -2.0
x_cur    = x_0
prob     = pos_tanh
p_cur    = prob(x_cur,-x_lim,x_lim)
print(qint.quad(prob,-x_lim,x_lim,args=(-x_lim,x_lim)))


x = np.arange(-3,3,0.01)
z = np.array([prob(xval,-x_lim,x_lim) for xval in x])

kept = []
for d,r in zip(Ds,Rs):
    x_trial = x_cur + d
    p_trial = prob(x_trial,-x_lim,x_lim)
    w       = p_trial/p_cur
    if w >= 1 or w > r:
        x_cur = x_trial
        n_accept += 1
        if n_accept > n_burn:
            kept.append(x_cur)

plt.plot(x,z,'k-',label='analytic')
plt.hist(kept,density=True,bins=50,color='#ff6666',label='MH realization')
plt.xlabel('x')
plt.ylabel('probability density')
plt.title('Metropolis-Hastings: $x_0$ = %s, $\delta$=%s, $f$=%2.2f'%(x_0,delta,n_accept/N_trials))
plt.legend()
plt.show()
print(n_accept)