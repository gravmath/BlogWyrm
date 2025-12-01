import cff_utils         as cff
import matplotlib.pyplot as plt
import numpy             as np
import pandas            as pd

#data from digitization of my image provided by https://apps.automeris.io/wpd/
nozzle_radii = np.array([
[3.63150867823765, 5.060975609756099],
[3.9252336448598126, 4.664634146341465],
[4.192256341789051, 4.176829268292685],
[4.432576769025367, 3.7500000000000018],
[4.65954606141522, 3.4146341463414647],
[4.939919893190922, 3.048780487804878],
[5.206942590120159, 2.8048780487804894],
[5.527369826435246, 2.560975609756099],
[5.887850467289718, 2.408536585365855],
[6.181575433911881, 2.3780487804878057],
[6.408544726301736, 2.3780487804879057],
[6.662216288384512, 2.408536585365855],
[6.969292389853138, 2.4695121951219523],
[7.222963951935915, 2.560975609756099],
[7.503337783711615, 2.6829268292682933],
[7.810413885180239, 2.835365853658537],
[8.09078771695594, 2.9573170731707332],
[8.411214953271026, 3.1402439024390247],
[8.691588785046727, 3.3079268292682933],
[8.98531375166889, 3.475609756097562],
[9.345794392523363, 3.7195121951219523],
[9.73297730307076, 3.9939024390243913],
[9.999999999999998, 4.146341463414635],
])

area = np.pi*nozzle_radii[:,1]**2

plt.plot(nozzle_radii[:,0],area)
plt.show()

Pe_P0 = 0.8
Astar = np.min(area)

#Quiz problem 5 Chapter 9
gamma   = 1.4
P0      = 400e3
T0      = 20 + 273
dthroat = 10
dexit   = 24 
A_Astar = (dexit/dthroat)**2
Mexit   = cff.inv_AMR(A_Astar,gamma,'supersonic')[0]
Pexit   = cff.P_P0(Mexit,gamma)*P0
print(A_Astar,Mexit,cff.P_P0(Mexit,gamma),Pexit)

#Quiz problem 8 Chapter 9
gamma   = 1.4
P0      = 400e3
T0      = 20 + 273
dthroat = 10
dexit   = 20 
A_Astar = (dexit/dthroat)**2
Mexit   = cff.inv_AMR(A_Astar,gamma,'supersonic')[0]
Pexit   = cff.P_P0(Mexit,gamma)*P0
Pdownst = cff.Pd_Pu(Mexit,gamma)*Pexit
print(A_Astar,Mexit,cff.P_P0(Mexit,gamma),Pexit,cff.Pd_Pu(Mexit,gamma),Pdownst)
