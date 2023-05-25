import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import levy_stable
from mpmath import *
import matplotlib as mpl
from matplotlib import rc
from datetime import datetime
import numba as nb
import h5py

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

#set font sizes
SMALL_SIZE = 22
MEDIUM_SIZE = 22

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=15)    # legend fontsize
################################################################################

f = h5py.File('scalar_data/scalar_data_s1/scalar_data_s1_p0.h5','r')
print(list(f.keys()))

dset  = f['tasks']
print(list(dset))
Ek    = dset['Ek'] 

#dset2 = f['scales']
#times = dset2['sim_time']
time = 0.01*np.arange(0,len(Ek[:,0,0]))
plt.plot(time,(Ek[:,0,0]),'o')
plt.plot(time,time,'k--')
#plt.plot(,'k--')

plt.figure(2)
dset  = f['tasks']
print(list(dset))
m     = dset['m']
plt.plot(time,(m[:,0,0]),'o')

plt.show()
