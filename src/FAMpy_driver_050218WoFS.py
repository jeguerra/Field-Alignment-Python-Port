#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 10:18:04 2022

Driver for FAMpy, loads some test data and performs alignment procedure

@author: jorge.guerra
"""

import numpy as np
import scipy as sp
import scipy.signal as sps
from matplotlib import cm
import matplotlib.pyplot as plt
import FAMpy_module as fpym
from netCDF4 import Dataset

if __name__ == '__main__':
       
       # User parameters
       tdex = 0
       zdex = 0
       fname = 'news-e_ENS_18_20180502_2100_2230'
       fiinput = '/media/jeguerra/BigDATA/work/verification_test_data/20180502/2100/' + fname + '.nc'
       
       # Field Alignment tuning
       niter = 1E4 # max iterations for spectral solver 
       lscale = 1 # change from 1~10, 1 is the best (according to...)
       vmode = 4 # viscosity or hyperviscosity model (2 or 4)
       pscale = 1.0 # conversion to pixel space for displacements      
       # Load test data and data dependent parameters
       m_fid = Dataset(fiinput, 'r', format='NETCDF4')
       
       # Pull in the 2D domain variables CHANGE AS NEEDED
       lons = m_fid.variables['xlon'][0,:]
       lats = m_fid.variables['xlat'][:,0]
       plons, plats = np.meshgrid(lons, lats)
       
       # Get the domain sizes
       nlon = lons.shape[0]
       nlat = lats.shape[0]
       
       # Pull the operating variable to align along with the number of members
       hmemAll =  m_fid.variables['comp_dz'][:] # reventador1
       
       #%% Rearrange data and apply windowing function QC
       hmem0 = hmemAll[0:5,:,:]
       osize = hmem0.shape 
       meml = range(osize[0])
       nmem = len(meml)
       hmemT = np.copy(hmem0)
       #'''
       win1 = sps.windows.tukey(nlon, alpha=0.25)
       win2 = sps.windows.tukey(nlat, alpha=0.25)
       win = np.outer(win1, win2)
       
       for imem in meml:
              hmemT_nominal = np.copy(hmemT[imem,:,:])
              hmemT[imem,:,:] = hmemT_nominal * win.T
              
       #%% Binarize fields? For object verification
       
       #%% INVOKE THE SOLVER!
       solver_args = (niter, lscale, vmode, nlon, nlat, nmem, meml, hmemT)
       hmemT_FA, hmem_FA, qxyT, qxy = fpym.computeAlignedEnsemble(*solver_args)
              
       #%% Make some nice plots
       
       fig, axs = plt.subplots(nrows=2, ncols=2)
       axs = axs.flatten()
       
       cl = 40.0 # dBZ for WoFS reflectivity
       xlims = (np.amin(lons), np.amax(lons)) # x window
       ylims = (np.amin(lats), np.amax(lats)) # y window
       prop_cycle = plt.rcParams['axes.prop_cycle']
       colors = prop_cycle.by_key()['color']
       cc = 0
       for pp in meml:
              
              if cc == len(colors):
                     cc = 0
                                          
              axs[0].contour(plons, plats, hmemT[pp,:,:], [cl], colors=colors[cc], linewidths=2.0)
              axs[0].set_xlim(*xlims)
              axs[0].set_ylim(*ylims)
              axs[1].contour(plons, plats, hmem_FA[:,:,pp].T, [cl], colors=colors[cc], linewidths=2.0)
              axs[1].set_xlim(*xlims)
              axs[1].set_ylim(*ylims)
              cc += 1
       
       AM_cf = axs[2].contourf(plons, plats, np.mean(hmemT[:,:,:], axis=0), 11, cmap=cm.jet)
       axs[2].set_xlim(*xlims)
       axs[2].set_ylim(*ylims)
       plt.colorbar(AM_cf, ax=axs[2])
       FAM_cf = axs[3].contourf(plons, plats, np.mean(hmem_FA[:,:,:], axis=2).T, 11, cmap=cm.jet)
       axs[3].set_xlim(*xlims)
       axs[3].set_ylim(*ylims)
       plt.colorbar(FAM_cf, ax=axs[3])
       plt.show()
       #plt.close()