import matplotlib        as mpl
import matplotlib.pyplot as plt
import numpy             as np
import pandas            as pd

from sklearn.cluster import KMeans

#generate the data
core      = pd.DataFrame()
core['x'] = np.random.normal(loc=0,scale=0.5,size=1000)
core['y'] = np.random.normal(loc=0,scale=1.0,size=1000)

fringe      = pd.DataFrame()
angs        = 2*np.pi*np.random.rand(1000)
rads        = 5 + 0.1*np.random.rand(1000)
fringe['x'] = rads*np.cos(angs)
fringe['y'] = rads*np.sin(angs)

tots        = pd.DataFrame()
tots['x']   = np.concatenate((core['x'],fringe['x']))
tots['y']   = np.concatenate((core['y'],fringe['y']))
#plt.plot(tots['x'],tots['y'],'k*',markersize=20)
plt.plot(core['x'],  core['y'],'r.')
plt.plot(fringe['x'],fringe['y'],'g.')
plt.show()

colors = ['red','green']
for q in range(100):
    kmeans       = KMeans(n_clusters=2, n_init=1).fit(tots)
    point_labels = kmeans.labels_
    fig_K  = plt.figure()
    ax_K   = fig_K.add_subplot(1,1,1)
    print(q,end=' ')
    for i in range(len(tots)):
        point = [tots.iloc[i][0],tots.iloc[i][1]]
        ax_K.plot(point[0],point[1],color=colors[point_labels[i]],marker='.',linewidth=0)
        ax_K.set_title('Iteration %s'%q)
    plt.show()
    del(fig_K)

colors = ['red','green','blue']
kmeans_cluster = KMeans(n_clusters=3, n_init=10).fit(tots)
point_labels = kmeans_cluster.labels_
fig_K  = plt.figure()
ax_K   = fig_K.add_subplot(1,1,1)
for i in range(len(tots)):
    point = [tots.iloc[i][0],tots.iloc[i][1]]
    ax_K.plot(point[0],point[1],color=colors[point_labels[i]],marker='.',linewidth=0)
plt.show()