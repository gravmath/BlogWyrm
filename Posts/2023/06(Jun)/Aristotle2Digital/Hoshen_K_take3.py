import matplotlib        as mpl
import matplotlib.pyplot as plt
import numpy             as np

#this implementation uses a modification of Gould & Tobochnik
def find_root(k,proper_labels):
    counter = 0
    while(proper_labels[k] != k and counter < len(proper_labels)):
        k = proper_labels[k]
        counter = counter + 1
    return k

def init_site_labels():
    proper_labels  = {-1:0,0:0}
    cluster_labels = {}
    for i in range(rows):
        for j in range(cols):
            cluster_labels[i*cols + j] = -1
    return cluster_labels, proper_labels  

def hoshen_kopleman(lattice):
    largest_label = 1
    rows, cols    = lattice.shape

    cluster_labels, proper_labels = init_site_labels()

    #forward sweep
    for i in range(rows):
        for j in range(cols):
            curr = lattice[i,j]
            if curr != 0:  #current cell occupied - only do the work if so
                up_index          = max(0,i-1)
                left_index        = max(0,j-1)
                up_linear_index   = up_index*cols + j     if i != 0 else 0
                left_linear_index = i*cols + left_index   if j != 0 else 0
                curr_linear_index = i*cols + j
                up                = lattice[up_index,j]   if i != 0 else 0
                left              = lattice[i,left_index] if j != 0 else 0
                if up == 0 and left == 0: #site is an island (so far)
                    proper_labels[largest_label]      = largest_label
                    cluster_labels[curr_linear_index] = largest_label
                    largest_label       += 1
                if up != 0 and left == 0: #site linked to the one above
                    cluster_labels[curr_linear_index] = cluster_labels[up_linear_index]
                if up == 0 and left != 0: #site linked to the one to the left
                    cluster_labels[curr_linear_index] = cluster_labels[left_linear_index]
                if up != 0 and left != 0: #site is a linker
                    cluster_label                     = min(cluster_labels[up_linear_index],cluster_labels[left_linear_index])
                    up_proper_label                   = proper_labels[cluster_labels[up_linear_index]]
                    left_proper_label                 = proper_labels[cluster_labels[left_linear_index]]
                    cluster_labels[curr_linear_index] = cluster_label
                    if up_proper_label > left_proper_label:
                        proper_labels[cluster_labels[up_linear_index]]   = proper_labels[cluster_labels[left_linear_index]]
                    else:
                        proper_labels[cluster_labels[left_linear_index]] = proper_labels[cluster_labels[up_linear_index]]

    #collapse the proper labels so that there is only one hop to the top of the tree
    collapsed_proper_labels = {}
    for k in proper_labels:
        collapsed_proper_labels[k] = find_root(k,proper_labels)

    return cluster_labels, collapsed_proper_labels

mode = "load"
rows = 4
cols = 7

M = np.array([[1,0,1,1,0,1,0],
              [1,1,0,0,1,1,1],
              [0,0,1,0,1,1,1],
              [1,1,1,1,1,1,1]])

rows = 7
cols = 7

M = np.array([[1,0,1,1,0,1,0],
              [1,1,0,0,1,1,1],
              [0,0,1,0,1,1,1],
              [1,1,1,1,1,1,1],
              [0,1,1,0,0,1,1],
              [1,0,0,0,1,1,1],
              [1,1,0,1,1,1,0]])

rows = 7
cols = 13
M = np.array([[0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0],
              [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1],
              [0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1],
              [0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0],
              [0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1],
              [1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1]])

cluster_labels, proper_labels = hoshen_kopleman(M)

#print(cluster_labels)
print(proper_labels)

set_unique_cluster_labels = set(proper_labels.values())
num_unique_clusters       = len(set_unique_cluster_labels)
translate = {}
for i,k in enumerate(set_unique_cluster_labels):
    translate[k] = i

Cs = np.zeros(rows*cols)

for i,k in enumerate(cluster_labels): 
    Cs[i] = translate[proper_labels[cluster_labels[k]]]

print(Cs.reshape(rows,cols))


cmp  = mpl.colors.ListedColormap(['w','#4472C4'])
fig = plt.figure()
ax  = fig.add_subplot(1,1,1)
ax.matshow(M,cmap=cmp)
ax.xaxis.set_major_locator(mpl.ticker.FixedLocator([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5]))
ax.yaxis.set_major_locator(mpl.ticker.FixedLocator([0.5,1.5,2.5,3.5,4.5,5.5,6.5]))
ax.grid('on')
plt.show()