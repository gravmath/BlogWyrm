import matplotlib        as mpl
import matplotlib.pyplot as plt
import matplotlib.dates  as md
import numpy             as np
import pandas            as pd

def calc_betaRP(factor,free_rate):
    return factor['beta']*(factor['P']-free_rate)

def calc_frac_change(arr):
    temp     = (np.roll(arr,-1) - arr)/arr
    temp     = np.roll(temp,1)
    temp[0]  = np.nan
    #temp[-1] = np.nan
    return temp
'''
x  = np.array([1,3,4,6,8,9,11,14])
y  = np.array([1,2,4,4,5,7,8,9])
n  = len(x)
Q  = np.polyfit(x,y,1)
Yf = Q[0]*x + Q[1]

fig = plt.figure()
ax  = fig.add_subplot(1,1,1)
ax.plot(x,Yf,'b-')

slopes = []
num_trials = 1000
for i in range(num_trials):
    x_noise = np.random.normal(0,0.1,n)
    y_noise = np.random.normal(0,0.5,n)
    Q = np.polyfit(x+x_noise,y+y_noise,1)
    slopes.append(Q[0])
    ax.plot(x+x_noise,y+y_noise,'r.')
ax.plot(x,y,'ko')
plt.show()

plt.hist(slopes)
plt.show()

#safety check
sx  = np.sum(x)
sx2 = np.sum(x*x)
sxy = np.sum(x*y)
sy  = np.sum(y)
sy2 = np.sum(y*y)
b = (sy*sx2 - sx*sxy)/(n*sx2 - sx*sx)
m = (n*sxy - sx*sy)/(n*sx2 - sx*sx)
'''

GDP = {'beta':0.6, 'P':0.07}
Inf = {'beta':0.8, 'P':0.05}
Gld = {'beta':-0.7,'P':0.08}
SP5 = {'beta':1.3, 'P':0.12}

factor_lst = [GDP,Inf,Gld,SP5]

free_rate  = 0.03
rate_asset = free_rate 
for factor in factor_lst:
    rate_asset += calc_betaRP(factor,free_rate)

print('rate asset value (%) = ',rate_asset*100.0)

df_NAS    = pd.read_csv('c:/Users/byecs/OneDrive/Documents/^IXIC.CSV')
df_MYT    = pd.read_csv('c:/Users/byecs/OneDrive/Documents/MMYT.CSV')
key_epochs= ['2012-01-01T00:00','2012-04-01T00:00','2012-07-01T00:00','2012-10-01T00:00',
             '2013-01-01T00:00','2013-04-01T00:00','2013-07-01T00:00','2013-10-01T00:00',
             '2014-01-01T00:00','2014-04-01T00:00','2014-07-01T00:00','2014-10-01T00:00']
key_mepochs = [md.date2num(np.datetime64(ke)) for ke in key_epochs]


fig_quotes = plt.figure(figsize=(8,4))
ax_quotes  = fig_quotes.add_subplot(1,1,1)
ax_quotes.plot(df_NAS['Date'],df_NAS['Adj Close']/np.mean(df_NAS['Adj Close']),label='NASDAQ')
ax_quotes.plot(df_MYT['Date'],df_MYT['Adj Close']/np.mean(df_MYT['Adj Close']),label='MakeMyTrip')
ax_quotes.xaxis.set_major_locator(mpl.dates.MonthLocator(interval=4))
ax_quotes.set_ylabel('Normalized Stock Price')
ax_quotes.legend()
plt.show()

df_NAS['frac'] = calc_frac_change(df_NAS['Adj Close'])
df_MYT['frac'] = calc_frac_change(df_MYT['Adj Close'])

plt.plot(df_NAS['Adj Close'],df_MYT['Adj Close'],color='#b1b198',marker='o',linewidth=0)
plt.xlabel('NASDAQ Adjusted Close ($)')
plt.ylabel('MyTripTime Adjusted Close ($)')
plt.show()

plt.plot(df_NAS['frac'],df_MYT['frac'],color='#b1b198',marker='o',linewidth=0)
plt.xlabel('NASDAQ Fractional Gain ($)')
plt.ylabel('MyTripTime Fractional Gain ($)')
plt.show()

print('Vaidya regression check: ',np.polyfit(df_NAS['frac'][1:],df_MYT['frac'][1:],1))
print('Mean of NASDAQ Composite',np.mean(df_NAS['Adj Close']))
print('Mean of MakeMyTrip',np.mean(df_MYT['Adj Close']))

#need to limit the range for the sample to avoid the NaN at the beginning
df_combined         = pd.DataFrame()
df_combined['Date'] = df_NAS['Date'][1:]
df_combined['NAS']  = df_NAS['frac'][1:]
df_combined['MYT']  = df_MYT['frac'][1:]

N_trials       = 10000
size_of_sample = 350
betas          = []
for n in range(N_trials):
    A = df_combined.sample(n=size_of_sample)
    Q = np.polyfit(A['NAS'],A['MYT'],1)
    betas.append(Q[0])

print(np.mean(betas),np.std(betas))

plt.plot(A['NAS'],A['MYT'],'ro')
plt.show()

fig_beta = plt.figure()
ax_beta  = fig_beta.add_subplot(1,1,1)
ax_beta.hist(betas,color='#a0a087')
ax_beta.set_xlabel('Value of $\\beta$')
ax_beta.set_ylabel('Relative Distribution')
ax_beta.yaxis.set_major_locator(mpl.ticker.FixedLocator(np.arange(0,4000,500)))
ax_beta.set_yticklabels(np.arange(0,4000,500)/N_trials)

plt.show()