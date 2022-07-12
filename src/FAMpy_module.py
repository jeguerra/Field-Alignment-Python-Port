#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 13:40:14 2021

Implementation of the field aligned mean (FAM) tendency
Module contains transforms, alignment solver, and advection

@author: jorge.guerra
"""

import sys
import cv2 as cv
import math as mt
import numpy as np
import scipy as sp
import scipy.sparse as sps
from matplotlib import cm
import matplotlib.pyplot as plt

def plotSurface(X, Y, Z):
       #fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
       #ax.plot_surface(X, Y, Z, cmap=cm.seismic,
       #         linewidth=0, antialiased=False)
       plt.figure(tight_layout=True)
       plt.contourf(X,Y,Z)
       plt.colorbar()
       plt.show()
       input('CHECKING SURFACE')       
       return

def highestPowerof2(n):
 
       res = 0
       for ii in range(n, 0, -1):
            
              # If i is a power of 2
              if ((ii & (ii - 1)) == 0):
               
                     res = ii
                     break
        
       return res

def numericalCleanUp(M):
       
       N = M.shape
       ZTOL = 1.0E-16
       MC = np.copy(M)
       # Clean up numerical zeros
       for ii in range(N[0]):
              for jj in range(N[1]):
                     if abs(M[ii,jj]) <= ZTOL:
                            MC[ii,jj] = 0.0
       return MC

def myImResize(IM, sz):
       
       szin = IM.shape
       
       xfrac = (szin[1] - 1) / sz[1]
       yfrac = (szin[0] - 1) / sz[0]
       
       pitchx = np.linspace(1+xfrac/2, szin[1]-xfrac/2, num=sz[1])
       pitchy = np.linspace(1+yfrac/2, szin[0]-yfrac/2, num=sz[0])
       
       x,y = np.meshgrid(pitchx,pitchy)
       
       bx = np.arange(1,szin[1]+1)
       by = np.arange(1,szin[0]+1)
       bdx,bdy = np.meshgrid(bx, by)
       
       out_spmf = sp.interpolate.RectBivariateSpline(bx, by, IM.T, kx=5, ky=5)
       out = out_spmf.ev(x.flatten(), y.flatten())
       out_2D = numericalCleanUp(np.reshape(out, x.shape, order='C'))
       
       #out = sp.interpolate.griddata((bdx.flatten(), bdy.flatten()), IM.flatten(), (x, y), method='cubic')
       #out = numericalCleanUp(out)
       
       return out_2D

def computeResizeDisplacementOutputs(X, Y, qtx, qty):
       
       # X input image
       # Y target image
       
       # Resize (by interpolation) to the next lower power of 2 on each dimension
       # This truncates, smooths, and forces uniform sampling
       
       sx = X.shape[0]
       sy = Y.shape[1]
       
       # Find the next lesser power of 2 to resize data arrays
       sx2 = highestPowerof2(sx)
       sy2 = highestPowerof2(sy)
       
       # Resample source and target images to the original size (sx, sy)
       qsavex = myImResize(qtx, (sx, sy))# * sy/sx2 
       qsavey = myImResize(qty, (sx, sy))# * sx/sy2 
       
       return qsavex, qsavey

# Advection by simple semi-Lagrangian transport (linear interpolation approach)
def advect(X0, qx, qy):
       
       if qx.size == 0:
              qxn = qx
              qyn = qy
              XOut_2D = np.copy(X0)
       else:
              qxn = np.nan_to_num(qx)
              qyn = np.nan_to_num(qy)
       
       bw = 50
       Xb = cv.copyMakeBorder(X0, bw, bw, bw, bw,cv.BORDER_REFLECT101)
       
       x0 = np.arange(1,X0.shape[1]+1)
       y0 = np.arange(1,X0.shape[0]+1)
       xx0, yy0 = np.meshgrid(x0, y0)
       
       xb = np.arange(1,Xb.shape[1]+1)
       yb = np.arange(1,Xb.shape[0]+1)
       xxb, yyb = np.meshgrid(xb, yb)
       
       xi = bw + xx0 - qxn
       yi = bw + yy0 - qyn
       
       # Apply interpolation to enforce displacements
       '''
       XOut_smpf = sp.interpolate.RectBivariateSpline(xb, yb, Xb.T, kx=1, ky=1)
       XOut = XOut_smpf.ev(xi.flatten(), yi.flatten())
       XOut_2D = numericalCleanUp(np.reshape(XOut, xx0.shape, order='C'))
       '''
       #'''
       # Linear interpolation is monotonic which is necessary for this advection scheme
       XOut_2D = sp.interpolate.griddata((xxb.flatten(), yyb.flatten()), Xb.flatten(),
                                      (xi, yi), method='linear', fill_value=0.0)
       #'''
       return XOut_2D

# Advection of an image by spatial derivatives
def advectOp(X0, qx, qy):
       
       if qx.size == 0:
              qxn = qx
              qyn = qy
              XOut = np.copy(X0)
       else:
              qxn = np.nan_to_num(qx)
              qyn = np.nan_to_num(qy)
              
              x0 = np.arange(1,X0.shape[1]+1)
              y0 = np.arange(1,X0.shape[0]+1)
              xx0, yy0 = np.meshgrid(x0, y0)
              
              # Advection by gradients (2nd order numpy gradient ops...)
              dXbx = np.gradient(X0, axis=1, edge_order=2)
              dXby = np.gradient(X0, axis=0, edge_order=2)
              dX = qxn * dXbx + qyn * dXby
              
              XOut = X0 - dX
       
       return XOut

def computeAlignedEnsemble(niter, lscale, vmode, nlon, nlat, nmem, meml, hmemT):
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
                     sx2 = highestPowerof2(sx)
                     sy2 = highestPowerof2(sy)
              
              XX.append(myImResize(thisHmemT, (sx2, sy2)))
       
       # Main solver loop
       qxyT_mem = np.zeros((nlat,nlon,2,nmem,nmem))
       
       for imem in meml:
              print('Reference member: ', str(imem))
              
              thisHmemT = hmemT[imem,:,:]
              
              # Compute displacement solver on upper triangle of data matrix
              for jmem in meml: #[imem+1:]:
                     if jmem != imem:
                            thatHmemT = hmemT[jmem,:,:]
                            
                            print('Aligning: ', str(imem), ' to ', str(jmem))
                            QXT, QYT = computeSolverSCA(XX[imem], XX[jmem], 
                                          niter, 0.1, lscale, vmode)
                            
                            qxTsave, qyTsave = computeResizeDisplacementOutputs(thisHmemT, thatHmemT, QXT, QYT) 
                            
                            qxyT_mem[:,:,0,imem,jmem] = qxTsave
                            qxyT_mem[:,:,1,imem,jmem] = qyTsave
              
              '''
              # Fill lower triangle by skew symmetry of displacements
              for jmem in meml[0:imem+1]:
                     qxyT_mem[:,:,0,imem,jmem] = -qxyT_mem[:,:,0,jmem,imem]
                     qxyT_mem[:,:,1,imem,jmem] = -qxyT_mem[:,:,0,jmem,imem]
              '''              
                            
              # Compute mean of all displacement vectors
              qxy = np.zeros((nlon,nlat,2))
              qxyT = np.zeros((nlat,nlon,2))
              qxyT[:,:,0] = np.mean(qxyT_mem[:,:,0,imem,:], axis=2)
              qxyT[:,:,1] = np.mean(qxyT_mem[:,:,1,imem,:], axis=2)
              qxy[:,:,0] = np.copy(qxyT[:,:,0].T)
              qxy[:,:,1] = np.copy(qxyT[:,:,1].T)
              
              print('Found displacement vector successfully.')
              
              # Now align the member imem by advection
              hmemT_FA = advect(thisHmemT, qxyT[:,:,0], qxyT[:,:,1])
              hmem_FA[:,:,imem] = hmemT_FA.T
              
       return hmemT_FA, hmem_FA, qxyT, qxy
       
def computeSolverSCA(X, Y, niter, wt, lscale, mode):
       
       # Yang and Ravela SCA method 2009 (Yang's Masters Thesis)
       # X input image
       # Y target image
       # niter number of total iterations
       # wt weird spectral truncation weight (= 1.0)
       # lscale part of weird spectral truncation coefficient
       # mode (= 2 or 4) for power law approximation type
       
       #%% User parameters and initializations
       convTol = 1.0E-4 # for displacements in pixel space
       errv0 = np.inf
       derrv = np.inf
       w1 = 1.0
       w2 = w1 / 3.0 # Newtonian fluid constants
       applyNormalization = True
       
       # Check that X and Y match size size
       if X.shape != Y.shape:
              print('ERROR: source and target image shapes do not match!')
              print('X, Y: ', X.shape, Y.shape)
              sys.exit()
       
       # Setup nominal size arrays
       szx = X.shape 
       nx = szx[0]
       ny = szx[1]
       pblmr1 = np.arange(1,nx+1)
       pblmr2 = np.arange(1,ny+1)
       
       # Compute the nominal pixel grid
       aa,bb = np.meshgrid(pblmr2,pblmr1)
       bdx, bdy = np.meshgrid(pblmr2,pblmr1)
       
       # Normalize the images prior to alignment - relative to source
       if applyNormalization:
           mxb = np.amin(X)
           dxb = np.amax(X) - mxb
           Xb = np.copy((X - mxb) / dxb)
           Y = np.copy((Y - mxb) / dxb)
       else:
           Xb = np.copy(X)

       # Set up the frequency domain
       fxdom = 2 * np.round(0.5 * (lscale * nx + nx))
       fydom = 2 * np.round(0.5 * (lscale * ny + ny))
       nzero = np.array([1.0E-8]) # near zero frequency...
       
       nfs1 = np.array(np.arange(-fxdom/2,0)) # negative frequencies
       pfs = np.array(np.arange(1,fxdom/2)) # positive frequencies
       ffx = np.concatenate((nfs1, nzero, pfs))
       
       nfs2 = np.array(np.arange(-fydom/2,0)) # negative frequencies
       pfs = np.array(np.arange(1,fydom/2)) # positive frequencies
       ffy = np.concatenate((nfs2, nzero, pfs))
       
       sdex1 = len(nfs1) # index to singularity x
       sdex2 = len(nfs2) # index to singularity y
       m, n = np.meshgrid(ffx, ffy)
       
       # Set up initial forcing and displacement arrays
       ux = np.zeros((nx,ny))
       uy = np.copy(ux)
       qtx = np.copy(ux)
       qty = np.copy(ux)
       
       #%% Main iteration loop
       ii = 1
       
       # Precompute deformation filters
       fac0 = w1 * (np.power(n, mode) + np.power(m, mode))
       fac1 = np.copy(w2 * np.power(n, 2.0) + fac0)
       fac2 = w2 * n * m
       fac4 = np.copy(w2 * np.power(m, 2.0) + fac0)
       
       xpy = fac1 - fac2
       ypx = fac4 - fac2
       zz1 = (fac1 * fac4 - fac2 * fac2) * (wt * np.reciprocal(fxdom))
       zz2 = (fac1 * fac4 - fac2 * fac2) * (wt * np.reciprocal(fydom))
       
       while (ii < niter) and (errv0 > convTol):
              
              # Compute difference p - q, interpolate and update
              qx = aa - qtx        
              qy = bb - qty
       
              #'''
              # Apply displacements by interpolation to updated coords (remapping)
              Xb_spmf = sp.interpolate.RectBivariateSpline(pblmr2, pblmr1, Xb.T, kx=1, ky=1)
              Xb_new = Xb_spmf.ev(qx.flatten(), qy.flatten())
              Xb = numericalCleanUp(np.reshape(Xb_new, qx.shape, order='C'))
              del(Xb_new)
              #'''
              '''
              Xb = sp.interpolate.griddata((bdx.flatten(), bdy.flatten()), \
                                           Xb.flatten(), (qx, qy), \
                                          method='linear', fill_value=0.0)
              '''
              '''
              Xb_new = advect(Xb, qx, qy)
              Xb = numericalCleanUp(Xb_new)
              del(Xb_new)
              '''                     
              # Compute the forcing
              dq = (Xb - Y)
              # Set all boundary edges to zero
              dq[0,:] = 0.0
              dq[:,0] = 0.0
              dq[-1,:] = 0.0
              dq[:,-1] = 0.0
              
              #print(ddx.shape, ddy.shape, X.shape)
              # Initiate force vectors and compute gradients (built in 2nd order numpy)
              dXbx = np.gradient(Xb, axis=1, edge_order=2)
              dXby = np.gradient(Xb, axis=0, edge_order=2)
              t1 = -dXbx * dq
              t2 = -dXby * dq
              
              # Compute spectral solution. Mode is 2 or 4
              # Yang and Ravela SCA method 2009 (Yang's Masters Thesis MIT)
              fFx = sp.fft.fftshift(sp.fft.fft2(t1, s=m.shape))
              fFy = sp.fft.fftshift(sp.fft.fft2(t2, s=n.shape))
              
              # Set values at the center frequency to avoid singularity...
              fFx[sdex2, sdex1] = 0.0
              fFy[sdex2, sdex1] = 0.0

              # Apply the deformation filters
              fux = (-fFx * xpy + fac2 * fFy) * np.reciprocal(zz1)
              fuy = (fFx * fac2 - ypx * fFy) * np.reciprocal(zz2)
              
              # Recover displacements
              ux = 0.5 * np.real(sp.fft.ifft2(sp.fft.ifftshift(fux), s=ux.shape))
              uy = 0.5 * np.real(sp.fft.ifft2(sp.fft.ifftshift(fuy), s=uy.shape))
                            
              # Change in the max-norm of displacements signals minimum
              errv1 = max(np.amax(np.abs(ux)), np.amax(np.abs(uy)))
              derrv = abs(errv1 - errv0)
              errv0 = np.copy(errv1)
              
              # update the solution
              qtx += ux
              qty += uy
              
              ii += 1
       
       qnorm = max(np.amax(np.abs(qtx)), np.amax(np.abs(qty)))
       print('Solution loop exit on: ', errv0, qnorm, ii)
       return qtx, qty