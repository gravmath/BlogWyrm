import matplotlib        as mpl
import matplotlib.pyplot as plt
import numpy             as np
import pandas            as pd

from sklearn.cluster import KMeans

#'read' height and weight data in
heights = np.array([65.0,73.0,59.0,61.0,75.0,67.0,68.0,70.0,62.0,66.0,77.0,75.0,74.0,70.0,61.0,58.0,66.0,59.0,68.0,61.0])
weights = np.array([220.0,160.0,110.0,120.0,150.0,240.0,230.0,220.0,130.0,210.0,190.0,180.0,170.0,210.0,110.0,100.0,230.0,120.0,210.0,130.0])

#determine the number of points
num_heights = len(heights)
num_weights = len(weights)
if num_weights == num_heights:
    num_pts = len(heights)
else:
    print('Data not equal in size')
    quit()

#normalize the data
n_heights = (heights - np.mean(heights))/np.std(heights)
n_weights = (weights - np.mean(weights))/np.std(weights)

#create a DataFrame
norm_data              = pd.DataFrame()
norm_data['n_heights'] = n_heights
norm_data['n_weights'] = n_weights

#make a list of points and point indices
K_means_centers = []
K_means_indices = []

#select first point at random and save the first element of index (index is list-like)
K_means_centers.append(norm_data.sample(1))
K_means_indices.append(K_means_centers[0].index)

#now loop over all points different from the current selected
dists = []
for i in range(num_pts):
    if i in K_means_indices:
        continue
    else:
        #calculate the distance
        dist = np.sqrt(norm_data.iloc[i][0]**2 + norm_data.iloc[i][1]**2)
        dists.append((i,dist))

colors = ['red','green','blue']
for q in range(100):
    kmeans       = KMeans(n_clusters=3, n_init=1).fit(norm_data)
    point_labels = kmeans.labels_
    fig_K  = plt.figure()
    ax_K   = fig_K.add_subplot(1,1,1)
    print(q,end=' ')
    for i in range(num_pts):
        point = [norm_data.iloc[i][0],norm_data.iloc[i][1]]
        ax_K.plot(point[0],point[1],color=colors[point_labels[i]],marker='o',linewidth=0)
        ax_K.set_title('Iteration %s'%q)
    plt.show()
    del(fig_K)