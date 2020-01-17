################################################################################################################################
# This python script reads in and computes differences between two simulations (AR and MO) on different temporal timescales: 
# 10, 30, 50, 100 and 200 years. Both simulations are 200 years long. The code is organized as follows:
# 1) First a function to compute mean and standard deviation is defined: MO_minus_AR.
# 2) Afterwards, a plotting function is defined:  Plot. This funtion will plot figures with a number of subplots depending 
#    on the timescale (i.e. Figure 1 will show 20 subplots each showing a 10-year mean, Figure 2 will show 6 subplots each 
#    showing a 30-year mean, and etc). The user can decide whether plotting using the polar or global projections.  
# 3) The different variables (SST, SIC, SAT, LW TOA, SW TOA, U, V and etc) are then uploaded - Users can choose which datset to
#    upload by commenting/uncommenting out the relevant bit of code
# 4) Finally the MO_minus_AR and Plot functions are called and the plots are produced. 
###############################################################################################################################
import iris 
import iris.plot as iplt
import iris.quickplot as qplt

import iris.analysis.cartography 
import cartopy.crs as ccrs
import cartopy.feature as cfe
from mpl_toolkits.basemap import Basemap, maskoceans, shiftgrid

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from iris.experimental.equalise_cubes import *
from iris.experimental.animate import*
import matplotlib.dates as mdates
import matplotlib.ticker as mtick

import pandas as pd
import scipy

########################## Function to compute mean and stdev of MO - AR differences ###########################################################################################

def MO_minus_AR(MO, AR):
 
 ##10 years differences 
 mean_10yr=iris.cube.CubeList()
 stdev_10yr=iris.cube.CubeList()
 stdev_MO_10yr=iris.cube.CubeList()

 for i in range (10, 200+1, 10):
   st=i-10
   mean_10yr.append( (MO[st:i,:,:] - AR[st:i,:,:]).collapsed('time', iris.analysis.MEAN) )
   stdev_10yr.append( (MO[st:i,:,:] - AR[st:i,:,:]).collapsed('time', iris.analysis.STD_DEV) )
   stdev_MO_10yr.append(MO[st:i,:,:].collapsed('time', iris.analysis.STD_DEV) )
  

 ##30 years differences 
 mean_30yr=iris.cube.CubeList()
 stdev_30yr=iris.cube.CubeList() 
 stdev_MO_30yr=iris.cube.CubeList()
 
 for i in range (30, 200+1, 30):
   st=i-30
   mean_30yr.append( (MO[st:i,:,:] - AR[st:i,:,:]).collapsed('time', iris.analysis.MEAN) )
   stdev_30yr.append( (MO[st:i,:,:] - AR[st:i,:,:]).collapsed('time', iris.analysis.STD_DEV) )
   stdev_MO_30yr.append(MO[st:i,:,:].collapsed('time', iris.analysis.STD_DEV) )


 ##50 years differences 
 mean_50yr=iris.cube.CubeList()
 stdev_50yr=iris.cube.CubeList()
 stdev_MO_50yr=iris.cube.CubeList()
 
 for i in range (50, 200+1, 50):
  st=i-50
  mean_50yr.append( (MO[st:i,:,:] - AR[st:i,:,:]).collapsed('time', iris.analysis.MEAN) )
  stdev_50yr.append( (MO[st:i,:,:] - AR[st:i,:,:]).collapsed('time', iris.analysis.STD_DEV) )
  stdev_MO_50yr.append(MO[st:i,:,:].collapsed('time', iris.analysis.STD_DEV) ) 


 ##100 years differences 
 mean_100yr=iris.cube.CubeList()
 stdev_100yr=iris.cube.CubeList()
 stdev_MO_100yr=iris.cube.CubeList()

 for i in range (100, 200+1, 100):
   st=i-100
   mean_100yr.append( (MO[st:i,:,:] - AR[st:i,:,:]).collapsed('time', iris.analysis.MEAN) )
   stdev_100yr.append( (MO[st:i,:,:] - AR[st:i,:,:]).collapsed('time', iris.analysis.STD_DEV) )
   stdev_MO_100yr.append(MO[st:i,:,:].collapsed('time', iris.analysis.STD_DEV) )
 

 ##200 years differences 
 mean_200yr= (MO[:,:,:] - AR[:,:,:]).collapsed('time', iris.analysis.MEAN) 
 stdev_200yr= (MO[:,:,:] - AR[:,:,:]).collapsed('time', iris.analysis.STD_DEV) 
 stdev_MO_200yr= MO[:,:,:].collapsed('time', iris.analysis.STD_DEV) 

 mean=iris.cube.CubeList([mean_10yr.merge_cube(), mean_30yr.merge_cube(), mean_50yr.merge_cube(), mean_100yr.merge_cube(), mean_200yr])
 stdev=iris.cube.CubeList([stdev_10yr.merge_cube(), stdev_30yr.merge_cube(), stdev_50yr.merge_cube(), stdev_100yr.merge_cube(), stdev_200yr])
 stdev_MO=iris.cube.CubeList([stdev_MO_10yr.merge_cube(), stdev_MO_30yr.merge_cube(), stdev_MO_50yr.merge_cube(), stdev_MO_100yr.merge_cube(), stdev_MO_200yr])

 
 return mean, stdev, stdev_MO


########################################### Plotting function ##########################################################################################################################

def Plot(mean, label, prj, lat0, long0, color_map, parallels, meridians, global_pr):

 for n in range (0, 4): 
 
  fig=plt.figure(tight_layout=True)
  if n == 0:  # set n_of_row and n_of_columns for subplots
   n_r = 5
   n_c = 4 
  if n == 1:  
   n_r = 2 
   n_c = 3
  if n == 2:  
   n_r = 2  
   n_c = 2 
  if n == 3:  
   n_r = 1  
   n_c = 2
  
  for i in range(0,mean[n].shape[0]):
  
   ax = fig.add_subplot(str(n_r), str(n_c), str(i+1)) #n_of_row, n_of_columns and plot_id 
   ax=plt.gca()

   if global_pr == 1:
    bmap=Basemap(projection= prj, llcrnrlat= -90,  urcrnrlat= 90, llcrnrlon=0,  urcrnrlon= 360, resolution='l')
   else: 
    bmap=Basemap(projection= prj,lat_0= lat0 ,lon_0= long0, resolution='l', ax=ax)
   
  
   lon= mean[n][i].coord('longitude').points
   lat= mean[n][i].coord('latitude').points
   x,y=bmap(*np.meshgrid(lon,lat))
    
   contours=bmap.contourf(x,y, mean[n][i].data, cmap=color_map ) 
   bmap.drawcoastlines()
   bmap.drawparallels(parallels)
   bmap.drawmeridians(meridians)
   cbar = plt.colorbar(contours, orientation='horizontal')

  if n == 0:  
   fig.suptitle( label + '  10-year mean', fontsize=16)
  if n == 1:  
   fig.suptitle( label + '  30-year mean', fontsize=16)
  if n == 2:  
   fig.suptitle(label + '  50-year mean', fontsize=16) 
  if n == 3:  
   fig.suptitle(label + '  100-year mean', fontsize=16)
 
 fig=plt.figure(figsize=(10,10))
 ax=plt.gca()

 if global_pr == 1:
   bmap=Basemap(projection= prj, llcrnrlat= -90,  urcrnrlat= 90, llcrnrlon=0,  urcrnrlon= 360, resolution='l')
 else: 
   bmap=Basemap(projection= prj,lat_0= lat0 ,lon_0= long0, resolution='l', ax=ax)
    
 lon= mean[4].coord('longitude').points
 lat= mean[4].coord('latitude').points
 x,y=bmap(*np.meshgrid(lon,lat))
 
 contours=bmap.contourf(x,y, mean[4].data, cmap=color_map) 
 bmap.drawcoastlines()
 bmap.drawparallels(parallels)
 bmap.drawmeridians(meridians)
 cbar = plt.colorbar(contours, orientation='horizontal') #format= '%.0e' #for P
 plt.title(label + '   200-year mean', fontsize=16)
 
 return

###################################################### LOADING DATA #################################################################################################################
#'''
#load annual SST data:
MO=iris.load_cube('/nerc/n02/n02/vittoria/seaice_annual_uar766_*.nc',  'sea surface temperature' ) 
AR_cubes=iris.load('/nerc/n02/n02/vittoria/seaice_annual_uas245_*.nc',     'sea surface temperature')
#concatenate cubes in one single dataset to plot
iris.util.unify_time_units(AR_cubes)
AR=AR_cubes.concatenate_cube()
#'''

'''
#load annual SIC data:
MO=iris.load_cube('/nerc/n02/n02/vittoria/seaice_annual_uar766_*.nc',  'ice area  (aggregate)' ) 
AR_cubes=iris.load('/nerc/n02/n02/vittoria/seaice_annual_uas245_*.nc',     'ice area  (aggregate)')
#concatenate cubes in one single dataset to plot
iris.util.unify_time_units(AR_cubes)
AR=AR_cubes.concatenate_cube()
AR.data=np.ma.masked_where(AR.data <= 0.15, AR.data)
MO.data=np.ma.masked_where(MO.data <= 0.15, MO.data)
'''

'''
#load Air temperature data 
MO=iris.load_cube('/nerc/n02/n02/vittoria/uar766_AirTemp_18502050.nc',  'air_temperature' ) 
AR=iris.load_cube('/nerc/n02/n02/vittoria/uas245_AirTemp_18502050.nc',     'air_temperature')
'''

'''
#load MSLP data 
MO=iris.load_cube('/nerc/n02/n02/vittoria/uar766_SLP_U_V_18502050.nc',  'air_pressure_at_sea_level' ) 
AR=iris.load_cube('/nerc/n02/n02/vittoria/uas245_SLP_U_V_18502050.nc',     'air_pressure_at_sea_level')
'''

'''
#load TOA LW flux data
MO=iris.load_cube('/nerc/n02/n02/vittoria/uar766_TOA_LW&SW_18502050.nc',  'toa_outgoing_longwave_flux' ) 
AR=iris.load_cube('/nerc/n02/n02/vittoria/uas245_TOA_LW&SW_18502050.nc',     'toa_outgoing_longwave_flux')
'''

'''
#load TOA SW flux data
MO=iris.load_cube('/nerc/n02/n02/vittoria/uar766_TOA_LW&SW_18502050.nc',  'toa_outgoing_shortwave_flux' ) 
AR=iris.load_cube('/nerc/n02/n02/vittoria/uas245_TOA_LW&SW_18502050.nc',     'toa_outgoing_shortwave_flux')
'''

'''
#load precipitation flux data
MO=iris.load_cube('/nerc/n02/n02/vittoria/uar766_PrecRate_18502050.nc',  'precipitation_flux' ) 
AR=iris.load_cube('/nerc/n02/n02/vittoria/uas245_PrecRate_18502050.nc',     'precipitation_flux')
'''

'''
#load U, V and compute wind speed 
U_MO=iris.load_cube('/nerc/n02/n02/vittoria/uar766_SLP_U_V_18502050.nc',  'x_wind' ) 
U_AR=iris.load_cube('/nerc/n02/n02/vittoria/uas245_SLP_U_V_18502050.nc',     'x_wind')
V_MO=iris.load_cube('/nerc/n02/n02/vittoria/uar766_SLP_U_V_18502050.nc',  'y_wind' ) 
V_AR=iris.load_cube('/nerc/n02/n02/vittoria/uas245_SLP_U_V_18502050.nc',     'y_wind')
MO= (U_MO**2. + V_MO**2.)**0.5
AR= (U_AR**2. + V_AR**2.)**0.5
'''

###################################### CALLING FUNCTIONS & PROCESS DATA ###############################################################################################################################

#use MO_minus_AR function to compute mean and stdev at each grid-points over different averaging periods 
mean, stdev = MO_minus_AR(MO, AR)

##calling plotting function:

##plot mean SST SH
Plot(mean, label='SST (DegC)', prj = 'ortho', lat0 = -90, long0 = 0, color_map='RdYlBu_r', parallels= np.arange(-80.,81.,20.), meridians = np.arange(-180.,181.,20.), global_pr=0 )
##plot mean SST NH
#Plot(mean, label='SST (DegC)',  prj = 'ortho', lat0 = 90, long0 = 0, color_map='RdYlBu_r', parallels= np.arange(-80.,81.,20.), meridians = np.arange(-180.,181.,20.), global_pr=0 )
##plot stdev SST SH
#Plot(stdev, label='SST',  prj = 'ortho', lat0 = -90, long0 = 0, color_map='RdYlBu_r', parallels= np.arange(-80.,81.,20.), meridians = np.arange(-180.,181.,20.), global_pr=0 )
##plot stdev SST NH
#Plot(stdev, label='SST',  prj = 'ortho', lat0 = 90, long0 = 0, color_map='RdYlBu_r', parallels= np.arange(-80.,81.,20.), meridians = np.arange(-180.,181.,20.), global_pr=0 )

##plot mean SIC SH
#Plot(mean, label='SIC', prj = 'ortho', lat0 = -90, long0 = 0, color_map='RdBu_r', parallels= np.arange(-80.,81.,20.), meridians = np.arange(-180.,181.,20.), global_pr=0)
##plot mean SIC NH
#Plot(mean, label='SIC', prj = 'ortho', lat0 = 90, long0 = 0, color_map='RdBu_r', parallels= np.arange(-80.,81.,20.), meridians = np.arange(-180.,181.,20.), global_pr=0)

##plot mean Air_Temp Global
#Plot(mean, label='SAT (DegC)', prj = 'gall', lat0 = None, long0 = None, color_map='RdYlBu_r', parallels= np.arange(-80.,81.,20.), meridians = np.arange(-180.,181.,20.), global_pr=1)

##plot mean TOA LW Global
#Plot(mean, label='LW TOA (W/m$^{2}$)', prj = 'gall', lat0 = None, long0 = None, color_map='YlOrRd', parallels= np.arange(-80.,81.,20.), meridians = np.arange(-180.,181.,20.), global_pr=1 )

##plot mean TOA SW Global
#Plot(mean, label='SW TOA (W/m$^{2}$)', prj = 'gall', lat0 = None, long0 = None, color_map='YlOrRd', parallels= np.arange(-80.,81.,20.), meridians = np.arange(-180.,181.,20.), global_pr=1  )

##plot mean Precipitation flux
#Plot(mean, label='P (kg/m$^{2}$/s)', prj = 'gall', lat0 = None, long0 = None, color_map='RdYlGn', parallels= np.arange(-80.,81.,20.), meridians = np.arange(-180.,181.,20.), global_pr=1 )

##plot mean wind speed at 10 m
#Plot(mean, label='10m Wind Speed', prj = 'gall', lat0 = None, long0 = None, color_map='bwr', parallels= np.arange(-80.,81.,20.), meridians = np.arange(-180.,181.,20.), global_pr=1 )
#Plot(mean, label='10m Wind Speed', prj = 'ortho', lat0 = 90, long0 = 0, color_map='bwr', parallels= np.arange(-80.,81.,20.), meridians = np.arange(-180.,181.,20.), global_pr=0 )

