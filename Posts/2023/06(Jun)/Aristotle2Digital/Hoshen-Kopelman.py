#Tamir Suliman implementation of the Hoshen-Kopelman algorithm for cluster labeling.
#The code implements a Union Find data structure which is 
#used to perform Hoshen-Kopelman algorithm on a binary matrix.
#
#https://tamirsuliman.medium.com/implementing-hoshen-kopelman-algorithm-using-python-30e144e8abf6


import numpy as np
import matplotlib.pyplot as plt
class UnionFind:
    def __init__(self, max_labels):
        self.labels = [0] * max_labels
        self.labels[0] = 0
        self.n_labels = max_labels
    def find(self, x):
        y = x
        while self.labels[y] != y:
            y = self.labels[y]
        while self.labels[x] != x:
            z = self.labels[x]
            self.labels[x] = y
            x = z
        return y
    def union(self, x, y):
        self.labels[self.find(x)] = self.find(y)
        return self.find(x)
    def make_set(self):
        self.labels[0] += 1
        assert self.labels[0] < self.n_labels
        self.labels[self.labels[0]] = self.labels[0]
        return self.labels[0]

def hoshen_kopelman(matrix, m, n):
    uf = UnionFind(m * n // 2)
    for i in range(m):
        for j in range(n):
            if matrix[i][j]:
                up = matrix[i - 1][j] if i > 0 else 0
                left = matrix[i][j - 1] if j > 0 else 0
                if up == 0 and left == 0:
                    matrix[i][j] = uf.make_set()
                elif up == 0 or left == 0:
                    matrix[i][j] = max(up, left)
                else:
                    matrix[i][j] = uf.union(up, left)
    return matrix

matrix1 = np.array([[1, 0, 0, 1, 0, 1],
                   [0, 1, 1, 0, 1, 1],
                   [1, 0, 0, 1, 1, 1],
                   [0, 1, 0, 0, 1, 1],
                   [0, 1, 1, 0, 1, 1],
                   [1, 0, 1, 0, 1, 1]])

print(hoshen_kopelman(matrix1, 6, 6))

hkm = hoshen_kopelman(matrix1, 6, 6)
fig, ax = plt.subplots()
min_val = 0
max_val = 6
m = 6 
n = 6
ax.matshow(hkm, cmap=plt.cm.Blues)
for i in range(m):
    for j in range(n):
        c = hkm[j, i]
        ax.text(i, j, str(c), va='center', ha='center')
plt.show()