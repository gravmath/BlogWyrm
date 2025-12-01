#scratch.py
import matplotlib        as mpl
import matplotlib.pyplot as plt
import numpy             as np
years = np.array([i for i in range(2021,2017,-1)])
'''
Possibly suspect data
#data from https://www.macrotrends.net/stocks/charts/SHEL/shell/gross-profit (and related)
BP_mkt_cap = 95250
EX_mkt_cap = 334369
CH_mkt_cap = 312465
SH_mkt_cap = 198403
BP = [44121,28207,49582,49319,38862,24627,23116,46438,63292,53216,68640,22858,55412]
EX = [31234,30942,53786,64695,53352,44475,53349,83480,89755,110658,106151,84805,63998]
CH = [66253,39705,62267,66894,53626,43483,56696,79759,81089,88767,88155,70270,54392]
SH = [72499,40354,70331,73847,59826,46917,43698,69804,72736,82656,85626,68482,54575]
'''
#From yahoo and confirmed via TDAmeritrade
#divide by a thousand to use $millions as the common value
BP_GP = np.array([22_860_000,10_184_000,27_583_000,28_880_000])/1000
CH_GP = np.array([48_309_000,24_475_000,30_534_000,44_905_000])/1000
EX_GP = np.array([64_886_000,8_127_000,55_958_000,67_733_000])/1000
SH_GP = np.array([35_849_000,-12_995_000,36_755_000,44_875_000])/1000
#net income common stockholder https://finance.yahoo.com/quote/XOM/financials?p=XOM
#include to be consistent with Politifact's '23 billion' value cited - https://www.politifact.com/factchecks/2022/mar/10/facebook-posts/yes-oil-companies-are-reporting-record-breaking-pr/
BP_NI = np.array([7_563_000,-20_306_000,4_025_000,9_382_000])/1000
CH_NI = np.array([15_625_000,-5_543_000,2_924_000,14_824_000])/1000
EX_NI = np.array([23_040_000,-22_440_000,14_340_000,20_840_000])/1000
SH_NI = np.array([20_101_000,-21_680_000,15_842_000,23_352_000])/1000
#as of close of market on 3/18/22
BP_mkt_cap = 94_174
EX_mkt_cap = 333_057
CH_mkt_cap = 314_977
SH_mkt_cap = 196_024
inflation  = np.array([1,1.018,1.031,1.079])

for y,i in zip(years,inflation):
    print(y,i)
    
bar_width  = 0.1

fig_GP    = plt.figure(figsize=(8,16))
ax_NI_GP  = fig_GP.add_subplot(2,1,1)
ax_NI_GP.xaxis.set_major_locator(mpl.ticker.FixedLocator(years))
ax_NI_GP.bar(years-0.15,BP_GP,label='BP',  width=bar_width,color='#f5f5dc',edgecolor='k',hatch='')
ax_NI_GP.bar(years-0.05,CH_GP,label='CVX', width=bar_width,color='#e4e4cb',edgecolor='k',hatch='/')
ax_NI_GP.bar(years+0.05,EX_GP,label='XOM', width=bar_width,color='#d3d3ba',edgecolor='k',hatch='\\')
ax_NI_GP.bar(years+0.15,SH_GP,label='SHEL',width=bar_width,color='#c2c2a9',edgecolor='k',hatch='/\\')
ax_NI_GP.legend(loc=9)
ax_NI_GP.yaxis.set_major_locator(mpl.ticker.FixedLocator([i*10_000 for i in range(-1,9)]))
ax_NI_GP.grid('on')
ax_NI_GP.set_ylabel('Gross Profit ($ million)')
ax_NI_GP.set_ylim([-15_000,80_000])

ax_I_GP  = fig_GP.add_subplot(2,1,2)
ax_I_GP.xaxis.set_major_locator(mpl.ticker.FixedLocator(years))
ax_I_GP.plot(years,inflation*BP_GP,label='BP   ($%s)'%BP_mkt_cap,color='#c2c2a9',marker='o',markeredgecolor='k',markersize=7)
ax_I_GP.plot(years,inflation*CH_GP,label='CVX  ($%s)'%CH_mkt_cap,color='#c2c2a9',marker='s',markeredgecolor='k',markersize=7)
ax_I_GP.plot(years,inflation*EX_GP,label='XOM  ($%s)'%EX_mkt_cap,color='#c2c2a9',marker='^',markeredgecolor='k',markersize=7)
ax_I_GP.plot(years,inflation*SH_GP,label='SHEL ($%s)'%SH_mkt_cap,color='#c2c2a9',marker='*',markeredgecolor='k',markersize=7)
ax_I_GP.legend()
ax_I_GP.yaxis.set_major_locator(mpl.ticker.FixedLocator([i*10_000 for i in range(-1,9)]))
ax_I_GP.grid('on')
ax_I_GP.set_ylabel('Inflation Adjusted Gross Profit ($ million)')
ax_I_GP.set_ylim([-15_000,80_000])
plt.show()

#https://www.politifact.com/factchecks/2022/mar/10/facebook-posts/yes-oil-companies-are-reporting-record-breaking-pr/
fig_NI    = plt.figure(figsize=(8,16))
ax_NI_NI  = fig_NI.add_subplot(2,1,1)
ax_NI_NI.xaxis.set_major_locator(mpl.ticker.FixedLocator(years))
ax_NI_NI.bar(years-0.15,BP_NI,label='BP',  width=bar_width,color='#f5f5dc',edgecolor='k',hatch='')
ax_NI_NI.bar(years-0.05,CH_NI,label='CVX', width=bar_width,color='#e4e4cb',edgecolor='k',hatch='/')
ax_NI_NI.bar(years+0.05,EX_NI,label='XOM', width=bar_width,color='#d3d3ba',edgecolor='k',hatch='\\')
ax_NI_NI.bar(years+0.15,SH_NI,label='SHEL',width=bar_width,color='#c2c2a9',edgecolor='k',hatch='/\\')
ax_NI_NI.legend(loc=9)
ax_NI_NI.yaxis.set_major_locator(mpl.ticker.FixedLocator([i*10_000 for i in range(-3,4)]))
ax_NI_NI.grid('on')
ax_NI_NI.set_ylabel('Net Income ($ million)')
ax_NI_NI.set_ylim([-30_000,30_000])

ax_I_NI  = fig_NI.add_subplot(2,1,2)
ax_I_NI.xaxis.set_major_locator(mpl.ticker.FixedLocator(years))
ax_I_NI.plot(years,inflation*BP_NI,label='BP   ($%s)'%BP_mkt_cap,color='#c2c2a9',marker='o',markeredgecolor='k',markersize=7)
ax_I_NI.plot(years,inflation*CH_NI,label='CVX  ($%s)'%CH_mkt_cap,color='#c2c2a9',marker='s',markeredgecolor='k',markersize=7)
ax_I_NI.plot(years,inflation*EX_NI,label='XOM  ($%s)'%EX_mkt_cap,color='#c2c2a9',marker='^',markeredgecolor='k',markersize=7)
ax_I_NI.plot(years,inflation*SH_NI,label='SHEL ($%s)'%SH_mkt_cap,color='#c2c2a9',marker='*',markeredgecolor='k',markersize=7)
ax_I_NI.legend()
ax_I_NI.yaxis.set_major_locator(mpl.ticker.FixedLocator([i*10_000 for i in range(-3,4)]))
ax_I_NI.grid('on')
ax_I_NI.set_ylabel('Inflation Adjusted Net Income ($ million)')
ax_I_NI.set_ylim([-30_000,30_000])
plt.show()
