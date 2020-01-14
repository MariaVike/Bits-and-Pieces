###############################################################################################################################
# This python script uses the IRIS package to upload and process simulation outputs stored as nc files. In the first part 
# of the code, the ENSO 3.4 index is computed starting from monthly SST data. In the second part, the ENSO index is plotted. 
#################################################################################################################################
import iris 
import iris.plot as iplt
import iris.quickplot as qplt

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from iris.experimental.equalise_cubes import *
from iris.experimental.animate import*

import pandas as pd
import scipy
import matplotlib.dates as mdates

'''
############################################## CALCULATE ENSO INDEX #############################################################
#upload data
MO_SST_cubes=iris.load('/nerc/n02/n02/vittoria/seaice_monthly_uar766_*.nc', iris.Constraint('sea surface temperature', latitude= lambda lat: -5 <= lat <= 5, longitude= lambda lon: -170 <= lon <= -120   ))

#concatenate cubes in one single dataset to plot
iris.util.unify_time_units(MO_SST_cubes)
MO_SST=MO_SST_cubes.concatenate_cube()

#compute area-weighted mean
MO_SST.coord('latitude').guess_bounds()
MO_SST.coord('longitude').guess_bounds()
area_weights_MO= iris.analysis.cartography.area_weights(MO_SST)

MO_SST_av = MO_SST.collapsed(['latitude', 'longitude'],
                                 iris.analysis.MEAN,
                                 weights=area_weights_MO)


#calculate the mean over a 30-year reference period
MO_SST_tmean=MO_SST_av[:30*12].collapsed('time', iris.analysis.MEAN)

#calculate anomalies:
ENSO_MO= MO_SST_av-MO_SST_tmean

#filter the time series with a 3-month running mean:
months_MO=ENSO_MO.shape[0]
ENSO_MO_3month=iris.cube.CubeList() 
for i in range (0, months_MO-3):
    ENSO_MO_3month.append(ENSO_MO[i:i+3].collapsed('time', iris.analysis.MEAN)) 
ENSO_MO_3month= ENSO_MO_3month.merge_cube()
print ENSO_MO_3month

iris.save(ENSO_MO_3month, 'ENSO_MO_18502050.nc')
'''

############################################## PLOT ENSO INDEX #######################################################################
#load data previously saved (for plotting purposes only, this way is faster)
ENSO_MO_3month=iris.load_cube('/nerc/n02/n02/vittoria/ENSO_MO_18502050.nc')

##select interval
#ENSO_MO_3month=ENSO_MO_3month[:50*12]

mean_MO=ENSO_MO_3month.collapsed('time', iris.analysis.MEAN)
stdev_MO=ENSO_MO_3month.collapsed('time', iris.analysis.STD_DEV)
print mean_MO.data, stdev_MO.data

months=ENSO_MO_3month.shape[0]
times= pd.date_range(start='1850-01-01', periods=months, freq='MS') #create array of dates with pandas 


fig = plt.figure()
ax = fig.add_subplot(211)

fig.autofmt_xdate()
plt.plot(times, ENSO_MO_3month.data, c='black',linewidth=1,  label='MONSOON')
plt.fill_between(times, 0.5, np.ma.masked_where(ENSO_MO_3month.data <= 0.5, ENSO_MO_3month.data) , alpha=0.5, facecolor='orangered')
plt.fill_between(times, -0.5, np.ma.masked_where(ENSO_MO_3month.data >= -0.5, ENSO_MO_3month.data) , alpha=0.5, facecolor='darkgreen')
plt.xlabel('Time(years)', fontsize=14)
plt.ylabel('ENSO Index', fontsize=14)
plt.title('ENSO Index MONSOON 1850-2050', fontsize=16)
plt.ticklabel_format(axis='y',useOffset=False) 
ax.xaxis.set_major_locator(mdates.YearLocator(10))
plt.legend()

plt.show()
