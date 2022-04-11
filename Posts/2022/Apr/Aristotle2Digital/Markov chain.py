#Markov chain.py
import matplotlib        as mpl
import matplotlib.pyplot as plt
import numpy             as np
import pandas            as pd

M = np.array([[0.9,0.5],[0.1,0.5]])
# N = M
# maxes = []
# mins  = []
# dets  = []
# print(M,np.linalg.det(M))
# for i in range(20):
#     N = np.linalg.matrix_power(M,i)
#     maxes.append(np.max(N))
#     mins.append(np.min(N))
#     dets.append(np.linalg.det(N))

# plt.semilogy(maxes)
# plt.show()

# plt.semilogy(mins)
# plt.show()

# plt.semilogy(dets)
# plt.show()
print(np.linalg.matrix_power(M,2))
print(np.linalg.matrix_power(M,10))
print(np.linalg.matrix_power(M,20))
print(np.linalg.det(np.linalg.matrix_power(M,20)))
#print(np.linalg.eigvals(M.T))

current_state = 0
states        = [current_state]
for i in range(20000):
    P = np.random.random()
    if current_state == 1:
        if P <= 0.9:
            current_state = 1
        else:
            current_state = 0
    if current_state == 0:
        if P <= 0.5:
            current_state = 1
        else:
            current_state = 0
    states.append(current_state)

days = pd.DataFrame()
days['weather'] = states

sample_sizes = [5,10,50,100,200,500,1000,2000,5000,10000,20000]

probs = []
for sample_size in sample_sizes:
    sample = np.array(days.sample(sample_size))
    sum    = np.sum(sample)
    print(sum)
    prob = sum/sample_size
    probs.append(prob)

print(sample_sizes,len(sample_sizes))
print(probs,len(probs))

fig = plt.figure()
ax  = fig.add_subplot(1,1,1)
ax.semilogx(sample_sizes,probs,'bo-')
ax.plot([1,20000],[0.833333,0.8333333],'k--')
ax.set_xlabel('Sample Size')
ax.set_ylabel('Probability of a Sunny Day')
ax.set_xlim([5,20000])
ax.set_ylim([0,1])
plt.show()

print(np.linalg.eig(M))

print(np.linalg.matrix_power(M,2))
print(np.linalg.matrix_power(M,10))
print(np.linalg.matrix_power(M,20))