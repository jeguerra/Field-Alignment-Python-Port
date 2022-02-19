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
       zdex = 9
       fname = 'reventador1'
       fiinput = '../test_data/' + fname + '.nc'
       
       # Field Alignment tuning
       niter = 1000 # max iterations for spectral solver 
       lscale = 1 #change from 1~10, 1 is the best (according to...)
       vmode = 4 #viscosity or hyperviscosity model (2 or 4)
       
       # Load test data and data dependent parameters
       m_fid = Dataset(fiinput, 'r', format='NETCDF4')
       
       # Pull in the 2D domain variables CHANGE AS NEEDED
       lons = m_fid.variables['longitude'][:]
       lats = m_fid.variables['latitude'][:]
       plons, plats = np.meshgrid(lons, lats)
       
       # Get the domain sizes
       nlon = lons.shape[0]
       nlat = lats.shape[0]
       
       # Pull the operating variable to align along with the number of members
       hmemAll =  m_fid.variables['Flight_levelA'][:] # reventador1
       #hmemAll =  m_fid.variables['__xarray_dataarray_variable__'][:,:,:,:] # kasatochi5
       
       # Perform alignment (double loop)
       
       # Rearrange data and apply windowing function QC
       hmem0 = hmemAll[-1,:,:,:]
       osize = hmem0.shape
       #meml = range(5) 
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
       #'''       
       # Initialize the aligned ensemble
       hmem_FA = np.zeros((nlon,nlat,nmem))
       
       XX = []
       # Resize images loop
       for imem in meml:
              thisHmemT = hmemT[imem,:,:]
              
              if imem == 0:
                     sx = thisHmemT.shape[0]
                     sy = thisHmemT.shape[1]
                     
                     # Find the next lesser power of 2
                     sx2 = fpym.highestPowerof2(sx)
                     sy2 = fpym.highestPowerof2(sy)
              
              XX.append(fpym.myImResize(thisHmemT, (sx2, sy2)))
       
       # Main solver loop 
       for imem in meml:
              print('Reference member: ', str(imem))
              qxyT_mem = np.zeros((nlat,nlon,2,nmem-1))
              
              thisHmemT = hmemT[imem,:,:]
              jj = 0
              for jmem in meml:
                     if jmem != imem:
                            thatHmemT = hmemT[jmem,:,:]
                            
                            print('Aligning: ', str(imem), ' to ', str(jmem))
                            QXT, QYT = fpym.computeSolverSCA(XX[imem], XX[jmem], 
                                          niter, 0.1, lscale, vmode)
                            
                            qxTsave, qyTsave = fpym.computeResizeDisplacementOutputs(thisHmemT, thatHmemT, QXT, QYT) 
                            
                            qxyT_mem[:,:,0,jj] = qxTsave
                            qxyT_mem[:,:,1,jj] = qyTsave
                            jj += 1
                            
              # Compute mean of all displacement vectors
              qxy = np.zeros((nlon,nlat,2))
              qxyT = np.zeros((nlat,nlon,2))
              qxyT[:,:,0] = np.mean(qxyT_mem[:,:,0,:], axis=2)
              qxyT[:,:,1] = np.mean(qxyT_mem[:,:,1,:], axis=2)
              qxy[:,:,0] = np.copy(qxyT[:,:,0].T)
              qxy[:,:,1] = np.copy(qxyT[:,:,1].T)
              
              print('Found displacement vector successfully.')
              
              # Now align the member imem by advection
              hmemT_FA = fpym.advect(thisHmemT, qxyT[:,:,0], qxyT[:,:,1])
              hmem_FA[:,:,imem] = hmemT_FA.T
              '''
              plt.contour(plons, plats, thisHmemT, [5.0E-4], colors='k', linewidths=2.0)
              plt.contour(plons, plats, hmemT_FA, [5.0E-4], colors='r', linewidths=2.0)
              plt.show()
              input()
              plt.close()
              '''
       #%% Make some nice plots
       
       fig, axs = plt.subplots(nrows=2, ncols=2)
       axs = axs.flatten()
       
       cl = 5.0E-4 # GOOD FOR THE REVENTADOR CASE ONLY!
       xlims = (-77.8, -77.4) # reventador1
       ylims = (-0.4, 0.0) # reventador1
       #xlims = (-176.0, -164.0) # kasatochi5
       #ylims = (46.5, 53.0) # kasatochi5
       prop_cycle = plt.rcParams['axes.prop_cycle']
       colors = prop_cycle.by_key()['color']
       cc = 0
       for pp in meml:
              
              if cc == len(colors):
                     cc = 0
                                          
              axs[0].contour(plons, plats, hmem0[pp,:,:], [cl], colors=colors[cc], linewidths=2.0)
              axs[0].set_xlim(*xlims)
              axs[0].set_ylim(*ylims)
              axs[1].contour(plons, plats, hmem_FA[:,:,pp].T, [cl], colors=colors[cc], linewidths=2.0)
              axs[1].set_xlim(*xlims)
              axs[1].set_ylim(*ylims)
              cc += 1
       
       AM_cf = axs[2].contourf(plons, plats, np.mean(hmem0[:,:,:], axis=0), 11, cmap=cm.jet)
       axs[2].set_xlim(*xlims)
       axs[2].set_ylim(*ylims)
       plt.colorbar(AM_cf, ax=axs[2])
       FAM_cf = axs[3].contourf(plons, plats, np.mean(hmem_FA[:,:,:], axis=2).T, 11, cmap=cm.jet)
       axs[3].set_xlim(*xlims)
       axs[3].set_ylim(*ylims)
       plt.colorbar(FAM_cf, ax=axs[3])
       plt.show()
       #plt.close()