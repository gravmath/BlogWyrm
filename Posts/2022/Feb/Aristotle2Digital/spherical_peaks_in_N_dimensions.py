import matplotlib        as mpl
import matplotlib.pyplot as plt
import numpy             as np
import scipy.special

def vol_coeff(n_dim):
    #https://www.usna.edu/Users/physics/mungan/_files/documents/Scholarship/HypersphereVolume.pdf
    return (np.pi)**(n_dim/2)/scipy.special.gamma(n_dim/2+1)

r          = 0.5
N_max      = 30
ds         = np.arange(1,N_max+1,1)
vol        = vol_coeff(ds)*r**ds
maj_lgtcks = [10**(-2*i) for i in range(0,15,1)]
min_lgtcks = [j*10**(-(i)) for i in range(0,15,1) for j in np.arange(0.1,1.0,0.1)]
#to show that log scale draws minor tick labels by default if they fall on base to some power
#min_lgtcks = []

fig = plt.figure()
ax  = fig.add_subplot(1,1,1)
ax.plot(ds,vol,color='#ff6666')
ax.tick_params(which='both',direction='in',top=True,right=True)
ax.set_yscale('log')
ax.set_xlim([1,20])
ax.set_ylim([1e-8,1])
ax.yaxis.set_minor_locator(mpl.ticker.FixedLocator(min_lgtcks))
ax.yaxis.set_major_locator(mpl.ticker.FixedLocator(maj_lgtcks))
#to show that log scale draws minor tick labels by default if they fall on base to some power
ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax.xaxis.set_major_locator(mpl.ticker.FixedLocator(ds[1::2]))
#to show that linear scale doesn't draw minor tick labels by default
ax.xaxis.set_minor_locator(mpl.ticker.FixedLocator(np.arange(1,N_max,0.2)))
#ax.xaxis.set_minor_locator(mpl.ticker.FixedLocator([])) 
ax.set_xlabel('Number of Dimensions ($N$)',fontsize=12)
ax.set_ylabel('Volume Fraction $V_{sphere}/V_{cube}$',fontsize=12)


plt.show()