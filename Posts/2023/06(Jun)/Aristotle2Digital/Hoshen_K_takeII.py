import matplotlib        as mpl
import matplotlib.pyplot as plt
import numpy             as np
import pdb 
import pickle            as pk
#https://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html
#a liar this man is as his 'This is a little bit sneaky.. ' section in his C 
#implementation is not discussed in the web page

#rather than use his implementation - this one uses a modification of Gould & Tobochnik
def init_site_labels(rows,cols):
    site_labels = {}
    for i in range(rows):
        for j in range(cols):
            site_labels[i*cols+j] = {'cluster':-1,'proper':-1,'site_index':[(i,j),i*cols+j]}
    return site_labels  

def hoshen_kopleman(lattice):
    largest_label = 1
    rows, cols    = lattice.shape
    site_labels   = init_site_labels(rows,cols)
    
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
                #print(i,j,curr,site_labels[i*cols+j],largest_label)
                if up == 0 and left == 0: #site is an island (so far)
                    #pdb.set_trace()
                    site_labels[curr_linear_index]['cluster'] = largest_label
                    site_labels[curr_linear_index]['proper']  = curr_linear_index
                    largest_label       += 1
                if up != 0 and left == 0: #site linked to the one above
                    site_labels[curr_linear_index]['cluster'] = site_labels[up_linear_index]['cluster']
                    site_labels[curr_linear_index]['proper']  = site_labels[up_linear_index]['proper']
                if up == 0 and left != 0: #site linked to the one to the left
                    site_labels[curr_linear_index]['cluster'] = site_labels[left_linear_index]['cluster']
                    site_labels[curr_linear_index]['proper']  = site_labels[left_linear_index]['proper']
                if up != 0 and left != 0: #site is a linker
                    #pdb.set_trace()
                    cluster_label = min(site_labels[up_linear_index]['cluster'],site_labels[left_linear_index]['cluster'])
                    proper_label  = min(site_labels[up_linear_index]['proper'], site_labels[left_linear_index]['proper'])
                    site_labels[curr_linear_index]['cluster'] = cluster_label
                    site_labels[curr_linear_index]['proper']  = proper_label
                    site_labels[up_linear_index]  ['proper']  = proper_label
                    site_labels[left_linear_index]['proper']  = proper_label
                #print(i,j,curr,site_labels[i*rows+j],largest_label)

    #backward sweep
    for i in range(rows-1,-1,-1):
        for j in range(cols-1,-1,-1):
            curr = lattice[i,j]
            if curr != 0:  #current cell occupied - only do work if so
                down_index          = min(rows-1,i+1)
                right_index         = min(cols-1,j+1)
                down_linear_index   = down_index*cols + j    if i != rows - 1 else 0
                right_linear_index  = i*cols + right_index   if j != cols - 1 else 0
                curr_linear_index   = i*cols + j
                down                = lattice[down_index,j]  if i != rows - 1 else 0
                right               = lattice[i,right_index] if j != cols - 1 else 0
                #pdb.set_trace()
                #print(i,j,curr,site_labels[i*cols+j],largest_label)
                if right == 0 and down == 0: #site is an island (so far)
                    continue
                if down != 0 and right == 0: #site linked to the one below
                    site_labels[curr_linear_index]['cluster'] = site_labels[down_linear_index]['cluster']
                    site_labels[curr_linear_index]['proper']  = site_labels[down_linear_index]['proper']
                if down == 0 and right != 0: #site linked to the one to the right
                    site_labels[curr_linear_index]['cluster'] = site_labels[right_linear_index]['cluster']
                    site_labels[curr_linear_index]['proper']  = site_labels[right_linear_index]['proper']
                if down != 0 and right != 0: #site is a linker
                    #pdb.set_trace()
                    cluster_label = min(site_labels[down_linear_index]['cluster'],site_labels[right_linear_index]['cluster'])
                    proper_label  = min(site_labels[down_linear_index]['proper'], site_labels[right_linear_index]['proper'])
                    site_labels[curr_linear_index]['cluster']  = cluster_label
                    site_labels[curr_linear_index]['proper']   = proper_label
                    site_labels[down_linear_index]  ['proper'] = proper_label
                    site_labels[right_linear_index]['proper']  = proper_label

    largest_label = np.max([site_labels[k]['cluster'] for k in site_labels])

    return site_labels, largest_label

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

site_labels, largest_label = hoshen_kopleman(M)

Cs = np.zeros(rows*cols)
for i,k in enumerate(site_labels): 
    print(i, site_labels[k])
    Cs[i] = site_labels[k]['cluster']

print(Cs.reshape(rows,cols))



cmp  = mpl.colors.ListedColormap(['w','b'])
fig = plt.figure()
ax  = fig.add_subplot(1,1,1)
ax.matshow(M,cmap=cmp)
ax.xaxis.set_major_locator(mpl.ticker.FixedLocator([0.5,1.5,2.5]))
ax.yaxis.set_major_locator(mpl.ticker.FixedLocator([0.5,1.5,2.5]))
ax.grid('on')
plt.show()