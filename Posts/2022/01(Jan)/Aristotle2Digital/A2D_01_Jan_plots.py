import matplotlib        as mpl
import matplotlib.pyplot as plt
import numpy             as np

x = np.arange(0,10.01,0.1)
f = np.exp(-x**2)
p = np.exp(-x)

fig = plt.figure()
ax  = fig.add_subplot(111)

ax.plot(x,f,'k-',label='$y(x) = e^{-x^2}$')
ax.plot(x,p,'r-',label='$y(x) = e^{-x}$')
ax.set_xlabel('$x$',fontsize=16)
ax.set_ylabel('$y(x)$',fontsize=16)
ax.legend(fontsize=12)
fig.savefig('D:/github/BlogWyrm/Posts/2022/Jan/Aristotle2Digital/importance_sampling.png')

plt.show()
