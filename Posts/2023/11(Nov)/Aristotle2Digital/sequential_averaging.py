import matplotlib.pyplot as plt
import numpy             as np
import scipy             as sp

def calc_seq_ave(curr_x,curr_ave,n):
    return (n-1)/n*curr_ave + 1/n*curr_x


N     = 1000
mu    = 11
sigma = 1.5
#samps = np.arang(N)

X = mu + sigma*np.random.randn(N)
curr_ave  = X[0]
curr_ave2 = X[0]**2
seq_aves  = [curr_ave]
seq_aves2 = [curr_ave2]
for i in range(1,N):
    curr_ave  = calc_seq_ave(X[i],curr_ave,i+1)
    curr_ave2 = calc_seq_ave(X[i]**2,curr_ave2,i+1)
    seq_aves.append(curr_ave)
    seq_aves2.append(curr_ave2)

fig  = plt.figure(figsize=(12,5))
axes = fig.subplots(1,2)
axes[0].plot(X,'b.')
axes[0].axhline(np.mean(X),color='k')
axes[0].set_xlabel('Sample')
axes[0].set_ylabel('X')
axes[1].plot(seq_aves,'ro-')
axes[1].axhline(np.mean(X),color='k')
axes[1].set_xlabel('Sample')
axes[1].set_ylabel('Sequential Average of X')
plt.show()

print(np.mean(X),np.std(X),curr_ave,np.sqrt(curr_ave2 - curr_ave**2))

Y = np.array([1,3,7,1,2,1,6])
curr_ave = Y[0]
seq_aves = [curr_ave]
for i in range(1,len(Y)):
    curr_ave = calc_seq_ave(Y[i],curr_ave,i+1)
    seq_aves.append(curr_ave)
#print(seq_aves)
