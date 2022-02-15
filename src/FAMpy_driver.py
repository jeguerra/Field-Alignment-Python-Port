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
import matplotlib.pyplot as plt
import FAMpy_module as fpym
from netCDF4 import Dataset

if __name__ == '__main__':
       
       # User parameters
       tdex = 0
       zdex = 9
       fname = 'reventador3'
       fiinput = '/scratch/HYSPLIT_DATA/reventador/' + fname + '.nc'
       
       # Field Alignment tuning
       niter = 10 # max iterations for spectral solver 
       lscale = 1 #change from 1~10, 1 is the best (according to...)
       vmode = 4 #hyperviscosity model
       
       meml = range(13)
       nmem = len(meml)
       
       # Load test data and data dependent parameters
       m_fid = Dataset(fiinput, 'r', format='NETCDF4')
       
       # Pull in the 2D domain variables CHANGE AS NEEDED
       lons = m_fid.variables['longitude'][:]
       lats = m_fid.variables['latitude'][:]
       
       # Get the domain sizes
       nlon = lons.shape[0]
       nlat = lats.shape[0]
       
       # Pull the operating variable to align
       hmemAll =  m_fid.variables['Flight_levelA'][:]
       osize = hmemAll.shape
       
       # Perform alignment (double loop)
       
       # Rearrange data and apply windowing function QC
       hmem0 = hmemAll[0,:,:,:]
       hmemT = hmem0.swapaxes(0,2)
       #hmemT = hmemT.swapaxes(0,1)
       #lons = lons.swapaxes(0,1)
       #lats = lats.swapaxes(0,1)
       
       win1 = sps.windows.tukey(nlon, alpha=0.25)
       win2 = sps.windows.tukey(nlat, alpha=0.25)
       win = np.outer(win1, win2)
       
       for imem in meml:
              hmemT_nominal = np.copy(hmemT[:,:,imem])
              hmemT[:,:,imem] = hmemT_nominal * win
              
       # Initialize the aligned ensemble
       hmem_FA = np.zeros((nlon,nlat,nmem))
       
       # Main solver loop 
       for imem in meml:
              print('Reference member: ', str(imem))
              qxyT_mem = np.zeros((nlat,nlon,2,nmem))
              
              for jmem in meml:
                     if jmem != imem:
                            print('Aligning: ', str(imem), ' to ', str(jmem))
                            qxyT_mem[:,:,0,imem], qxyT_mem[:,:,1,imem] = \
                                   fpym.computeSolverSCA(hmemT[:,:,imem], hmemT[:,:,jmem], 
                                                         niter, 0.1, lscale, vmode)
                                   
              # Compute mean of all displacement vectors
              qxy = np.zeros((nlon,nlat,2))
              qxyT = np.zeros((nlon,nlat,2))
              qxyT[:,:,0] = np.mean(qxyT_mem[:,:,0,:], axis=3)
              qxyT[:,:,1] = np.mean(qxyT_mem[:,:,1,:], axis=3)
              qxy[:,:,0] = np.copy(qxyT[:,:,0].T)
              qxy[:,:,1] = np.copy(qxyT[:,:,1].T)
              
              print('Found displacement vector successfully.')
              
              # Now align the member imem by advection
              hmemT_FA = fpym.advect(hmemT[:,:,imem], qxyT[:,:,0], qxyT[:,:,1])
              hmem_FA[:,:,imem] = hmemT_FA.T
              
       # Make some nice plots