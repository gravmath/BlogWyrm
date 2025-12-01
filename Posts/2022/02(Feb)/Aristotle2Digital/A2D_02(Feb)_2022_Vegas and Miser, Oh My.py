import matplotlib       as mpl
import matplotlib.pyplot as plt
import numpy             as np
import scipy.special     as spec

def vol_nsphere(ndim,R):
    coeff      = np.pi**(ndim/2)/spec.gamma(ndim/2+1)
    vol_sphere = coeff*R**(ndim)
    return vol_sphere

L         = 1
R         = 0.5
Ns        = np.arange(1,21,1)
vol_boxes = 1**(Ns)
vol_fracs = vol_nsphere(Ns,R)/vol_boxes

fig = plt.figure(figsize=(5,5))
fig.set_facecolor('#ffeeee')

#create and tailor the axes
ax  = fig.add_subplot(1,1,1)
ax.plot(Ns,vol_fracs,color='#ff6666',linestyle='-',linewidth=2)

#tailor the x-axis
ax.xaxis.set_tick_params(which='both',direction='in',top=True,right=True)
ax.xaxis.set_major_locator(mpl.ticker.FixedLocator(np.arange(1,21,2)))
ax.xaxis.set_minor_locator(mpl.ticker.FixedLocator(np.arange(1,21.02,0.2)))
ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax.set_xlim([1,20])
ax.set_xlabel('Number of Dimesions $N_D$')

#tailor the y-axis
ax.set_yscale('log')
ax.yaxis.set_tick_params(which='both',direction='in',top=True,right=True)
ax.yaxis.set_major_locator(mpl.ticker.FixedLocator(np.array([(0.01)**(i) for i in range(0,10,1)])))
y_min_tcks = np.array([j*10**(i) for i in range(-8,0,1) for j in range(1,10,1)])
ax.yaxis.set_minor_locator(mpl.ticker.FixedLocator(y_min_tcks))
ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax.set_ylim([1e-8,1])
ax.set_ylabel('Volume Fraction')

print(vol_fracs[-1])
#plt.show()
#fig.savefig('vol_frac.png',dpi=857/5)