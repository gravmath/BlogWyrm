import matplotlib           as mpl
import matplotlib.pyplot    as plt
import mpl_toolkits.mplot3d as mpl3
import matplotlib.cm        as cm
import matplotlib.ticker    as mpt
import numpy                as np

def d1(x,y,z,mu):
    return np.sqrt( (x+mu)**2 + y**2 + z**2 )

def d2(x,y,z,mu):
    return np.sqrt( (x-1+mu)**2 + y**2 + z**2)

def U_Jacobi(x,y,z,mu,cutoff):
    u = 0.5*(x**2+y**2) + (1-mu)/d1(x,y,z,mu) + mu/d2(x,y,z,mu)
    u[np.where(u>cutoff)] = cutoff
    return u

sun_mass   = 1.989e30
earth_mass = 5.972e24
moon_mass  = 7.350e22

print( (earth_mass+moon_mass)/(sun_mass+earth_mass+moon_mass))
print( moon_mass/(earth_mass+moon_mass))

mu      = 0.5 #Feb UTH & Mar UTH
mu      = 0.05 #Mar UTH
#mu = 0.01
xmin, xmax, ymin, ymax = -1.5,1.5,-1.5,1.5
delta_x = 0.001
delta_y = 0.001
z_val   = 0.0
x = np.arange(xmin,xmax+delta_x,delta_x)
y = np.arange(ymin,ymax+delta_y,delta_y)
X, Y = np.meshgrid(x,y)
Z    = U_Jacobi(X,Y,z_val,mu,4.5)
'''
# Create a figure and a 3D axis
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# Plot the surface
Uxy = ax.plot_surface(X, Y, Z, cmap='Greens_r',edgecolor='none', linewidth=0, antialiased=False)
ax.set_zlim([-1,5])
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_zlabel(r'$U(x,y,%s)$'%z_val)
# Add a color bar
fig.colorbar(Uxy)#
# Show the plot
plt.show()
'''

CJs_half = [3.0,2.5, 2.0,1.75,
            1.72,1.71,1.70,1.65,
            1.6,1.526,1.4863,1.45]  #C_J/2 = U

num_rows = 3
num_cols = 4
fig, axes  = plt.subplots(nrows=num_rows,ncols=num_cols,figsize=(8,8))
axes = axes.flatten()
for i, ax in enumerate(axes):
    #con = ax.contour(X,Y,Z,levels=[1.375,1.4,1.5,1.6,1.7,1.8,1.9])
    ax.set_aspect('equal', adjustable='box')
    ax.tick_params(
        axis   = 'both',     # Apply to both x and y axes
        which  = 'both',     # Apply to major and minor ticks
        bottom = True,       # Remove ticks on bottom
        top    = True,       # Remove ticks on top
        left   = True,       # Remove ticks on left
        right  = True,       # Remove ticks on right
        labelbottom = False, # Remove x-axis labels
        labelleft   = False  # Remove y-axis labels
    )
    ax.set_facecolor('#f0f0f0')
    con = ax.contourf(X,Y,Z,levels=[CJs_half[i],100])
    cur_col = i % num_cols
    cur_row = (i-cur_col)/num_cols % num_rows
    if cur_row == 2:
        ax.set_xlabel(r'$x$')
        ax.tick_params(
        labelbottom = True, # Include x-axis labels
        )
    if cur_col == 0:
        ax.set_ylabel(r'$y$')
        ax.tick_params(
        labelleft   = True, # Include y-axis labels
        )

    ax.set_title(r'$C_J$ = '+'%s'% (2*CJs_half[i]) )
#fig.colorbar(con)
plt.tight_layout()
plt.show()