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

import computeDerivativeMatrix as derv

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
       
       out_spmf = sp.interpolate.RectBivariateSpline(bx, by, IM.T)
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
       
       # Find the next lesser power of 2
       sx2 = highestPowerof2(sx)
       sy2 = highestPowerof2(sy)
       
       # Resample source and target images to the original size (sx, sy)
       qsavex = myImResize(qtx, (sx, sy)) * sy/sx2 
       qsavey = myImResize(qtx, (sx, sy)) * sx/sy2 
       
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
       
       bw = 25
       Xb = cv.copyMakeBorder(X0, bw, bw, bw, bw,cv.BORDER_REFLECT101)
       
       x0 = np.arange(1,X0.shape[1]+1)
       y0 = np.arange(1,X0.shape[0]+1)
       xx0, yy0 = np.meshgrid(x0, y0)
       
       xb = np.arange(1,Xb.shape[1]+1)
       yb = np.arange(1,Xb.shape[0]+1)
       xxb, yyb = np.meshgrid(xb, yb)
       
       xi = bw + xx0 + qxn
       yi = bw + yy0 + qyn
       
       # Apply interpolation to enforce displacements
       '''
       XOut_smpf = sp.interpolate.RectBivariateSpline(xb, yb, Xb.T)
       XOut = XOut_smpf.ev(xi.flatten(), \
                           yi.flatten())
       XOut_2D = numericalCleanUp(np.reshape(XOut, xx0.shape, order='C'))
       '''
       # Linear interpolation is monotonic which is necessary for this advection scheme
       XOut_2D = sp.interpolate.griddata((xxb.flatten(), yyb.flatten()), Xb.flatten(),
                                      (xi, yi), method='linear')
       
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
       
       # Compute CS derivative matrices on the index grid
       '''
       ddx, d2dx = derv.computeCubicSplineDerivativeMatrix(x0, False, False, False, False, None)
       ddy, d2dy = derv.computeCubicSplineDerivativeMatrix(y0, False, False, False, False, None) 
       dX = qxn * (ddx @ X0.T).T + qyn * (ddy @ X0)
       '''
       #'''
       dXbx = np.gradient(X0, axis=1, edge_order=2)
       dXby = np.gradient(X0, axis=0, edge_order=2)
       dX = qxn * dXbx + qyn * dXby
       #'''
       XOut = X0 + dX
       
       return XOut
       
def computeSolverSCA(X, Y, niter, wt, lscale, mode):
       
       # Yang and Ravela SCA method 2009 (Yang's Masters Thesis)
       # X input image
       # Y target image
       # niter number of total iterations
       # wt weird spectral truncation weight (= 1.0)
       # lscale part of weird spectral truncation coefficient
       # mode (= 2 or 4) for power law approximation type
       
       #%% User parameters and initializations
       convTol = 1.0E-6
       errv0 = np.inf
       derrv = np.inf
       w1 = 1.0
       w2 = w1 / 3.0 # Newtonian fluid constants
       
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
       mxb = np.amin(X)
       dxb = np.amax(X) - mxb
       Xb = np.copy((X - mxb) / dxb)
       Y = np.copy((Y - mxb) / dxb)

       # Set up the frequency domain
       lscale = 2 * np.round(0.5 * (lscale * nx + nx))
       nzero = np.array([1.0E-7]) # near zero frequency...
       nfs = np.array(np.arange(-lscale/2,0)) # negative frequencies
       pfs = np.array(np.arange(1,lscale/2)) # positive frequencies
       ffs = np.concatenate((nfs, nzero, pfs))
       lffs = len(ffs)
       m, n = np.meshgrid(ffs, ffs)
       
       # Set up initial forcing and displacement arrays
       ux = np.zeros((nx,ny))
       uy = np.copy(ux)
       qtx = np.copy(ux)
       qty = np.copy(ux)
       
       #%% Main iteration loop
       ii = 1
       
       # Compute CS derivative matrices on the index grid
       ddx, d2dx = derv.computeCubicSplineDerivativeMatrix(pblmr2, False, False, False, False, None)
       ddy, d2dy = derv.computeCubicSplineDerivativeMatrix(pblmr1, False, False, False, False, None)
       
       #ddx = derv.computeCompactFiniteDiffDerivativeMatrix1(pblmr2, 4)
       #ddy = derv.computeCompactFiniteDiffDerivativeMatrix1(pblmr1, 4)
       
       ddxs = sps.csr_matrix(ddx); del(ddx)
       ddys = sps.csr_matrix(ddy); del(ddy)
       
       pdex = np.ix_(pblmr1-1,pblmr2-1)
       while (ii < niter) and (derrv > convTol):
              
              # Index local displacement to regions
              ssx = ux[pdex]
              ssy = uy[pdex]
              
              # Compute difference p - q, interpolate and update
              qx = aa - ssx        
              qy = bb - ssy
              qtx += ssx 
              qty += ssy

              # Apply displacements by interpolation to updated coords (remapping)
              Xb_spmf = sp.interpolate.RectBivariateSpline(pblmr2, pblmr1, Xb.T)
              Xb_new = Xb_spmf.ev(qx.flatten(), qy.flatten())
              Xb = numericalCleanUp(np.reshape(Xb_new, qx.shape, order='C'))
              del(Xb_new)
              
              #Xb = sp.interpolate.griddata((bdx.flatten(), bdy.flatten()), \
              #                             Xb.flatten(), (qx, qy), \
              #                            method='linear', fill_value=0.0)
                                   
              # Compute the forcing
              dq = (Xb - Y)
              # Set all boundary edges to zero
              dq[0,:] = 0.0
              dq[:,0] = 0.0
              dq[-1,:] = 0.0
              dq[:,-1] = 0.0
              
              #print(ddx.shape, ddy.shape, X.shape)
              # Initiate force vectors and compute gradients
              #dXbx = np.gradient(Xb, axis=1, edge_order=2)
              #dXby = np.gradient(Xb, axis=0, edge_order=2)
              dXbx = (ddxs.dot(Xb.T)).T
              dXby = ddys.dot(Xb)
              t1 = -dXbx * dq
              t2 = -dXby * dq
                            
              # Compute spectral solution. Mode is 2 or 4
              # Yang and Ravela SCA method 2009 (Yang's Masters Thesis MIT)
              fFx = sp.fft.fftshift(sp.fft.fft2(t1, s=(lffs, lffs)))
              fFy = sp.fft.fftshift(sp.fft.fft2(t2, s=(lffs, lffs)))
              
              # Set values at the center frequency to avoid singularity...
              lscalei = len(nfs)
              fFx[lscalei, lscalei] = 0.0
              fFy[lscalei, lscalei] = 0.0

              # Apply the deformation filters
              fac0 = w1 * (np.power(n, mode) + np.power(m, mode))
              fac1 = np.copy(w2 * np.power(n, 2.0) + fac0)
              fac2 = w2 * n * m
              fac4 = np.copy(w2 * np.power(m, 2.0) + fac0)
              
              xpy = fac1 - fac2
              ypx = fac4 - fac2
              zz = (fac1 * fac4 - fac2 * fac2) * (wt / lscale)
                            
              fux = (-fFx * xpy + fac2 * fFy) * np.reciprocal(zz)
              fuy = (fFx * fac2 - ypx * fFy) * np.reciprocal(zz)
              
              # Recover displacements
              ux = np.real(sp.fft.ifft2(sp.fft.ifftshift(fux)))
              uy = np.real(sp.fft.ifft2(sp.fft.ifftshift(fuy)))
              
              #plotSurface(aa, bb, ux[pdex])
              
              # Change in the max-norm of displacements signals minimum
              errv1 = max(np.amax(np.abs(ux)), np.amax(np.abs(uy)))
              derrv = abs(errv1 - errv0)
              errv0 = np.copy(errv1)
              #print(np.amax(np.abs(ux)), np.amax(np.abs(uy)))
              
              ii += 1
              
       print('Solution loop exit on: ', derrv, ii)
       return qtx, qty