import matplotlib        as mpl
import matplotlib.pyplot as plt
import numpy             as np
import scipy             as sp

Da    = 1.66053906660e-27 # kg
k_B   = 1.380649e-23      #J⋅K^−1
N_A   = 6.02214076e23     # per mol

def MB_dist(v,m,T):
    therm_E = k_B*T
    K_E     = m*v**2/2.0
    coeff   = 4*np.pi*(m/(2*np.pi*therm_E))**(1.5)
    expon   = -K_E/therm_E
    return coeff*v**2*np.exp(expon)

def calc_Vrms(T,M):
    return np.sqrt(3*k_B*T/M)

def calc_Vmp(T,M):
    return np.sqrt(2*k_B*T/M)    

def calc_Vave(T,M):
    return np.sqrt(8*k_B*T/M/np.pi)    

M_N2  = 2*14.007*Da       #mass of nitrogen in kg
M_air = 28.97             #g/mol
T     = 293
v_mp  = calc_Vmp(T,M_N2) 
f_mp  = MB_dist(v_mp,M_N2,T)
v_rms = calc_Vrms(T,M_N2)
f_rms = MB_dist(v_rms,M_N2,T)
v_ave = calc_Vave(T,M_N2)
f_ave = MB_dist(v_ave,M_N2,T)

print(v_mp,v_ave,v_rms)
print(np.sqrt(2),np.sqrt(8/np.pi),np.sqrt(3))
print(v_ave/v_mp,v_rms/v_mp)

V = np.arange(0,2000,10)
F = MB_dist(V,M_N2,T)

plt.plot(V,F,'k-')
plt.xlabel('Molecular Speed $v$ (km/s)',fontsize=16)
plt.ylabel('$f(v)$',fontsize=16)
plt.ylim([0,0.0023])
plt.vlines(v_mp, 0,MB_dist(v_mp, M_N2,T),color='blue', label='$v_{mp}$')
plt.vlines(v_ave,0,MB_dist(v_ave,M_N2,T),color='cyan', label='$v_{ave}$')
plt.vlines(v_rms,0,MB_dist(v_rms,M_N2,T),color='red',  label='$v_{rms}$')
plt.legend(fontsize=16)
plt.show()

Ts = np.array([50,100,200,300,500,1000])
Cs = ['indigo','blue','cyan','green','orange','red']
for T,C in zip(Ts,Cs):
    F = MB_dist(V,M_N2,T)
    plt.plot(V,F,label='T = %s'%T,color=C)
plt.xlabel('Molecular Speed $v$ (km/s)',fontsize=16)
plt.ylabel('$f(v)$',fontsize=16)
plt.ylim([0,0.006])
plt.legend(fontsize=16)
plt.show()    