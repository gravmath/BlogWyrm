import matplotlib        as mpl
import matplotlib.pyplot as plt
import numpy             as np
import pandas            as pd


#Matthew Stafford stats from Wikipedia
years  = [2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021]
GP     = [10,3,16,16,16,16,16,16,16,16,8,16,17]
yards  = [2267,535,5038,4967,4650,4257,4262,4327,4446,3777,2499,4084,4886]
TDs    = [13,6,41,20,29,22,32,24,29,21,19,26,41]
Prtg   = [61,91.3,97.2,79.8,84.2,85.7,97,93.3,99.3,89.9,106,96.3,102.9]
Place  = [29,np.nan,5,22,19,21,9,13,6,25,6,14,6]
wins   = [2,6,10,4,7,11,7,9,9,6,3,5,12]
losses = [14,10,6,12,9,5,9,7,7,10,12,11,5]
ties   = [0,0,0,0,0,0,0,0,0,0,1,0,0]

df = pd.DataFrame()
df['Years'] = years
df['GP']    = GP
df['yards'] = yards
df['TDs']   = TDs
df['Prtg']  = Prtg
df['Place'] = Place
df['Ws']    = wins
df['Ls']    = losses
df['Ts']    = ties

fig  = plt.figure(figsize=(10,5))
ax0  = fig.add_subplot(2,2,1)
ax0.plot(df['Years'][df['GP']>10],df['Prtg'][df['GP']>10],'bo') 
ax1  = fig.add_subplot(2,2,2)
ax1.plot(df['Prtg'][df['GP']>10],df['Ws'][df['GP']>10],'bo')
ax2  = fig.add_subplot(2,2,3)
ax2.plot(df['Years'][df['GP']>10],df['Place'][df['GP']>10],'yo')
ax3  = fig.add_subplot(2,2,4)
ax3.plot(df['Place'][df['GP']>10],df['Ws'][df['GP']>10],'yo')
print(df['Years'][df['GP']>10])
plt.show()