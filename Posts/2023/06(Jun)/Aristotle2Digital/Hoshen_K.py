import copy
import matplotlib        as mpl
import matplotlib.pyplot as plt
import numpy             as np
import pdb 
import pickle            as pk
#https://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html
#a liar this man is as his 'This is a little bit sneaky.. ' section in his C 
#implementation is not discussed in the web page

def init_labels(rows,cols):
    labels = np.arange(0,rows*cols)
    return labels  

def find_root(x,labels):
    counter = 0
    while(labels[x] != x and counter < len(labels)):
        x = labels[x]
        counter = counter + 1
    return x

def find(x,labels):
    root = find_root(x,labels)
    
    parent = labels[x]
    counter = 0
    while(labels[x] != x and counter < len(labels)):
        labels[x] = root
        x         = parent
        counter   = counter + 1
    return root

def union(x,y,labels):
    if x > y:
        labels[find(x,labels)] = labels[find(y,labels)]
    else:
        labels[find(y,labels)] = labels[find(x,labels)]

def hoshen_kopleman(lattice):
    largest_label = 1
    rows, cols    = lattice.shape
    labels        = init_labels(rows,cols)
    
    for i in range(rows):
        for j in range(cols):
            curr = lattice[i,j]
            if curr != 0:  #current cell occupied - only do the work if so
                linear_index      = i*cols + j
                up_index          = max(0,i-1)
                left_index        = max(0,j-1)
                up_linear_index   = up_index*cols + j   if i != 0 else 0
                left_linear_index = i*cols + left_index if j != 0 else 0
                up                = lattice[up_index,j]   if i != 0 else 0
                left              = lattice[i,left_index] if j != 0 else 0
                pdb.set_trace()
                print(i,j,curr,linear_index,up_index,left_index,up_linear_index,left_linear_index,up,left)
                if up == 0 and left == 0: #site is an island (so far)
                    labels[linear_index] = largest_label
                    largest_label       += 1
                if up != 0 and left == 0: #site linked to the one above
                    labels[linear_index] = find(up_linear_index,labels)
                if up == 0 and left != 0: #site linked to the one to the left
                    labels[linear_index] = find(left_linear_index,labels)                    
                if up != 0 and left != 0: #site is a bridge
                    union(up_linear_index,left_linear_index,labels)
                    labels[linear_index] = find(up_linear_index,labels)
                print(i,j,curr,linear_index,up_index,left_index,up_linear_index,left_linear_index,up,left)                    
                print('labels = ',labels)
                print('************************\n\n')


    return labels, largest_label

mode = "load"
rows = 4
cols = 4
cmp  = mpl.colors.ListedColormap(['k','w'])
if mode == "create":
    M    = np.random.randint(0,2,size=(rows,cols))
    with open('C:/Users/byecs/Downloads/temp/lattice.pkl', 'wb') as lattice_file:
        pk.dump(M, lattice_file)
elif mode == "load":
    with open('C:/Users/byecs/Downloads/temp/lattice7.pkl','rb') as lattice_file:
        M = pk.load(lattice_file)
        N = copy.deepcopy(M)

labels, largest_label = hoshen_kopleman(M)

Ls = labels.reshape(rows,cols)
print(M,'\n*******\n',labels,'\n*******\n',largest_label,'\n*******\n',Ls*M)

fig = plt.figure()
ax  = fig.add_subplot(1,1,1)
ax.matshow(N,cmap=cmp)
ax.xaxis.set_major_locator(mpl.ticker.FixedLocator([0.5,1.5,2.5]))
ax.yaxis.set_major_locator(mpl.ticker.FixedLocator([0.5,1.5,2.5]))
ax.grid('on')
plt.show()