# Code to recreate Caesar 2021 plot--------

# Import modules that are used
import numpy as np # mathematical functions
import pandas as pd
import scipy
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt # use Matplotlib to plot the data
from netCDF4 import Dataset
import os
# Doing the Lowess filtering
import statsmodels.api as sm # install using pip install statsmodels
lowess = sm.nonparametric.lowess   
from scipy.interpolate import griddata
import xlrd
# Writting the outputs as csv files
import csv
#from matplotlib.image import pil_to_array
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
from PIL import Image

os.environ['PROJ_LIB'] = r'/users/research/mupton/3. RSL_North_Atlantic/Caesar2021paper_plots'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

#%% Working directory
os.chdir(r'/users/research/mupton/3. RSL_North_Atlantic/Caesar2021paper_plots')

#%% load data from Spoon et al., submitted - %T quinqueloba data to argue for changes in SPG productivity 
wb = xlrd.open_workbook('LC_outputs/Spooner et al.xls') 
sheet = wb.sheet_by_index(0) #1st sheet:Data Set S5: RAPID-17-5P faunal count data. Core collected at 61.4817o N, 19.5360o W, 2303 m depth.
Tquin_time = sheet.col_values(0,3)[0:69] # smoothed mean SS (/mu m)
Tquin = sheet.col_values(2,3)[0:69]
Tquin_low = sheet.col_values(4,3)[0:69]
Tquin_high = sheet.col_values(3,3)[0:69]

Tquin2_time = sheet.col_values(6,3) # smoothed mean SS (/mu m)
Tquin2 = sheet.col_values(8,3)
Tquin2_low = sheet.col_values(10,3)
Tquin2_high = sheet.col_values(9,3)


Tquin.extend(Tquin2)
Tquin_low.extend(Tquin2_low[0:45])
Tquin_high.extend(Tquin2_high[0:45])
Tquin_time.extend(Tquin2_time[0:45])
Tquin_time = np.asarray(Tquin_time, dtype=float)
Tquin = np.asarray(Tquin, dtype=float)
Tquin_low = np.asarray(Tquin_low, dtype=float)
Tquin_high = np.asarray(Tquin_high, dtype=float)

Tquin_time_lin = np.linspace(int(Tquin_time[0]),int(Tquin_time[-1]),int(Tquin_time[0]-Tquin_time[-1]+1))
Tquin_lin = griddata(Tquin_time,Tquin, Tquin_time_lin, method='cubic')
Tquin_lowess_50y = lowess(Tquin_lin, Tquin_time_lin, 0.031, it=5 ) # 100/1621*50=3.1%
Tquin_low_lin = griddata(Tquin_time, Tquin_low, Tquin_time_lin, method='cubic')
Tquin_low_lowess_50y = lowess(Tquin_low_lin, Tquin_time_lin, 0.031, it=5 ) # 100/1621*50=3.1%
Tquin_high_lin = griddata(Tquin_time, Tquin_high, Tquin_time_lin, method='cubic')
Tquin_high_lowess_50y = lowess(Tquin_high_lin, Tquin_time_lin, 0.031, it=5 ) # 100/1621*50=3.1%

mean_spoon = pd.DataFrame(Tquin_lowess_50y, columns = ['Year','AMOC'])
lwr_spoon = pd.DataFrame(Tquin_low_lowess_50y, columns = ['Year','AMOC_lwr'])
upr_spoon = pd.DataFrame(Tquin_high_lowess_50y, columns = ['Year','AMOC_upr'])

spoon_df_lwr = pd.merge(mean_spoon,lwr_spoon)
spoon_df = pd.merge(spoon_df_lwr,upr_spoon)
#Spooner output dataframe
spoon_df.to_csv('LC_outputs/lowess_outputs/Spooner.csv')


#%% load data from Thornalley et al.  2018, Fig. 3e
wb = xlrd.open_workbook('LC_outputs/Thornalley_Figure_3.xls') 
sheet = wb.sheet_by_index(4) #1st sheet: AMOC indixces
# 48 sediment core
jpc48 = sheet.col_values(9,2) # smoothed mean SS (/mu m)
jpc48_time = sheet.col_values(8,2)

jpc48 = np.asarray(jpc48[0:68], dtype=float)
jpc48_time = np.asarray(jpc48_time[0:68], dtype=float)
jpc48_time_lin = np.linspace(int(jpc48_time[0]),int(jpc48_time[-1]),int(jpc48_time[0]-jpc48_time[-1]+1))
jpc48_lin = griddata(jpc48_time, jpc48, jpc48_time_lin, method='cubic')
#jpc48_lowess_20y = lowess(jpc48_lin, jpc48_time_lin, 0.012, it=5) # 100/1615*20 = 1.2%
jpc48_lowess_50y = lowess(jpc48_lin, jpc48_time_lin,0.031, it=5) # 100/1615*50 = 3.1%
# # thorn48 output dataframe
mean_thorn48_df = pd.DataFrame(jpc48_lowess_50y, columns = ['Year','AMOC'])
mean_thorn48_df.to_csv('LC_outputs/lowess_outputs/thorn48.csv')

# 56 sediment core
jpc56 = sheet.col_values(4,3) # smoothed mean SS (/mu m)
jpc56_time = sheet.col_values(2,3)

jpc56 = np.asarray(jpc56[0:92], dtype=float)
jpc56_time = np.asarray(jpc56_time[0:92], dtype=float)

jpc56_time_lin = np.linspace(int(jpc56_time[0]),int(jpc56_time[-1]),int(jpc56_time[0]-jpc56_time[-1]+1))
jpc56_lin = griddata(jpc56_time, jpc56, jpc56_time_lin, method='cubic')
#jpc56_lowess_20y = lowess(jpc56_lin, jpc56_time_lin, 0.038, it=5) # 100/530*20 = 3.8%
jpc56_lowess_50y = lowess(jpc56_lin, jpc56_time_lin,0.098 , it=5) # 100/530*50 = 9.8%
# # thorn56 output dataframe
mean_thorn56_df = pd.DataFrame(jpc56_lowess_50y, columns = ['Year','AMOC'])
mean_thorn56_df.to_csv('LC_outputs/lowess_outputs/thorn56.csv')



#%% load data from Sherwood et al. 2011, Fig. 5
wb = xlrd.open_workbook('LC_outputs/Sherwood et al Fig 5 data.xls') 
sheet = wb.sheet_by_index(0) 

d15N_time_ave = sheet.col_values(3,3)
d15N_time_min = sheet.col_values(1,3)
d15N_time_max = sheet.col_values(2,3)
d15N = sheet.col_values(6,3) # Bulk delta 15N (pro mill)
d15N_std = sheet.col_values(7,3)  # Bulk delta 15N std (pro mill)

d15N_time_ave = np.asarray(d15N_time_ave[0:10], dtype=float)
d15N_time_min = np.asarray(d15N_time_min[0:10], dtype=float)
d15N_time_max = np.asarray(d15N_time_max[0:10], dtype=float)
d15N = np.asarray(d15N[0:10], dtype=float) # Bulk delta 15N (pro mill)
d15N_std = np.asarray(d15N_std[0:10], dtype=float)  # Bulk delta 15N std (pro mill) 

d15N2_time = np.asarray(sheet.col_values(0,16), dtype=float)
d15N2 = np.asarray(sheet.col_values(1,16), dtype=float) # Bulk delta 15N (pro mill)
d15N2_low = np.asarray(sheet.col_values(2,16), dtype=float)  # Bulk delta 15N std (pro mill)
d15N2_high = np.asarray(sheet.col_values(3,16), dtype=float)

d15N2_lowess_20y = lowess(d15N2, d15N2_time, 0.26, it=5) # 100/77*20 = 26%
d15N2_low_lowess_20y = lowess(d15N2_low, d15N2_time, 0.26, it=5) # 100/77*20 = 26%
d15N2_high_lowess_20y = lowess(d15N2_high, d15N2_time, 0.26, it=5) # 100/77*20 = 26%

mean_sher = pd.DataFrame(d15N2_lowess_20y, columns = ['Year','AMOC'])
lwr_sher = pd.DataFrame(d15N2_low_lowess_20y, columns = ['Year','AMOC_lwr'])
upr_sher = pd.DataFrame(d15N2_high_lowess_20y, columns = ['Year','AMOC_upr'])

sher_df_lwr = pd.merge(mean_sher,lwr_sher)
sher_df = pd.merge(sher_df_lwr,upr_sher)
#Sherwood output dataframe
sher_df.to_csv('LC_outputs/lowess_outputs/Sherwood2011.csv')


#%% load data from Thibodeau et al. 2018, Fig. 5--------
wb = xlrd.open_workbook('LC_outputs/Thibodeau et al data.xls') 

sheet = wb.sheet_by_index(1)

MD99_time = sheet.col_values(0,1)
MD99 = sheet.col_values(1,1) # delta 18O (pro mill)
MD99_time = np.asarray(MD99_time[0:179], dtype=float)
MD99 = np.asarray(MD99[0:179], dtype=float) # delta 18O (pro mill)
MD99_time_lin = np.linspace(int(MD99_time[0]),int(MD99_time[-1]),int(MD99_time[0]-MD99_time[-1]+1))
MD99_lin = griddata(MD99_time, MD99, MD99_time_lin, method='cubic')
MD99_lowess_20y = lowess(MD99_lin, MD99_time_lin, 0.04, it=5) # 100/1253*50 = 4%
# Thibodeau 2018 md99 output dataframe
mean_thibmd99_df = pd.DataFrame(MD99_lowess_20y, columns = ['Year','AMOC'])
mean_thibmd99_df.to_csv('LC_outputs/lowess_outputs/thibmd99_2018.csv')

#%% load temperature-based AMOC index by Rahmstorf et al., 2015 and Caesar et al., 2018
mann = np.loadtxt('LC_outputs/AMOC_Mann.txt')
mann_lowess_50y = lowess(mann[:,1], mann[:,0], 0.046, it=5) # 100%/1096*50 = 4.56%
mann_upper = mann[:,1] + mann[:,2]
mann_lowess_upper_50y = lowess(mann_upper, mann[:,0], 0.046, it=5) # 100%/1096*50 = 4.56%

mann_lower = mann[:,1] - mann[:,2]
mann_lowess_lower_50y = lowess(mann_lower, mann[:,0], 0.046, it=5) # 100%/1096*50 = 4.56%

#Ramhsoft outputdatafraome
mean_rahm = pd.DataFrame(mann_lowess_50y, columns = ['Year','AMOC'])
lwr_rahm = pd.DataFrame(mann_lowess_lower_50y, columns = ['Year','AMOC_lwr'])
upr_rahm = pd.DataFrame(mann_lowess_upper_50y, columns = ['Year','AMOC_upr'])

rahm_df_lwr = pd.merge(mean_rahm,lwr_rahm)
rahm_df = pd.merge(rahm_df_lwr,upr_rahm)

rahm_df.to_csv('LC_outputs/lowess_outputs/rahm2015.csv')


plt.show()

