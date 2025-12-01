import matplotlib        as mpl
import matplotlib.pyplot
import numpy             as np
import scipy.special     as spec

plt  = mpl.pyplot
tick = mpl.ticker

#define a function to compute the 
#volume of a hypersphere, of radius
# $r$ in n-dimensions $ndim$
def sphere_vol(radius=0.5,ndim=3):
    coeff = np.pi**(ndim/2)/spec.gamma(ndim/2+1)
    return coeff*radius**ndim

#make some data for figure 1
X1        = np.arange(0,5.1,0.1)
Y1        = 2**X1
X1tcks    = np.arange(0,5,1)
X1mintcks = np.arange(0,5.5,0.5)
Y1tcks    = np.array([2**i for i in range(6)])
Y1mintcks = np.arange(2,36,2)

#now make figure 1
  #create the figure and add the axes
fig1 = plt.figure()
ax1  = fig1.add_subplot(1,1,1)
  #plot the data
ax1.plot(X1,Y1,color='#ff6666',linestyle='--')
  #setup the look-and-feel of the ticks
ax1.axes.tick_params(axis='both',which='both',direction='in',\
                     bottom=True,top=True,left=True,right=True)
ax1.axes.tick_params(axis='both',which='major',grid_color='k',\
                     grid_linestyle='--')
ax1.axes.tick_params(axis='y',which='minor',grid_color='gray',\
                     grid_linewidth=0.5,grid_linestyle=':')
  #tailor the y-axis
ax1.set_yscale('log')
ax1.yaxis.set_major_locator(tick.FixedLocator(Y1tcks))
ax1.yaxis.set_major_formatter(tick.FixedFormatter(Y1tcks))
ax1.yaxis.set_minor_locator(tick.FixedLocator(Y1mintcks))
ax1.yaxis.set_minor_formatter(tick.NullFormatter())
  #set the state of the grid and limit the axes scope
ax1.axes.grid(which='both')
ax1.set_xlim([0,5])
ax1.set_ylim([1,32])

plt.show()

#make some data common to figures 2 and 3
ndim_max       = 21
X23            = np.arange(1,ndim_max,1)
Y23            = sphere_vol(ndim=X23)
Y23_tcks       = [10**(-2*i) for i in range(0,int(np.ceil(ndim_max/2)),1)]
Y23_labels     = [str(i) for i in Y23_tcks]
Y23_min_tcks   = [j*10**(-i) for i in range(0,int(np.ceil(ndim_max/2)),1) \
                  for j in np.arange(0.9,0.0,-0.1)]

#now make figure 2
  #create the figure and add the axes
fig2 = plt.figure()
ax2  = fig2.add_subplot(1,1,1)
  #plot the data
ax2.plot(X23,Y23,color='#ff6666')
  #setup the look-and-feel of the ticks
ax2.tick_params(which='minor',direction='in',top=True,right=True)
ax2.tick_params(which='major',direction='inout',top=True,right=True)
  #tailor the x-axis
ax2.xaxis.set_major_locator(tick.FixedLocator(X23))
every_nth = 2
for n, label in enumerate(ax2.xaxis.get_ticklabels()):
    if (n+1) % every_nth != 0:
        label.set_visible(False)
ax2.xaxis.set_minor_locator(tick.FixedLocator(np.arange(1,ndim_max,0.2)))
ax2.xaxis.set_minor_formatter(tick.NullFormatter())
ax2.set_xlabel('Number of Dimensions ($N$)',fontsize=12)
  #tailor the y-axis
ax2.set_yscale('log')
ax2.yaxis.set_major_locator(tick.FixedLocator(Y23_tcks))
ax2.yaxis.set_minor_locator(tick.FixedLocator(Y23_min_tcks))
ax2.yaxis.set_minor_formatter(tick.NullFormatter())
ax2.set_ylabel('Volume Fraction $V_{sphere}/V_{cube}$',fontsize=12)
  #turn on the grid and limit the axes scope
ax2.axes.grid(which='major')
ax2.set_xlim([1,ndim_max])
ax2.set_ylim([1e-8,1])
plt.show()