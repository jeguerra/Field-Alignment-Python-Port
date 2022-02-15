#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 08:59:22 2019

@author: TempestGuerra
"""

import scipy.linalg as scl
import numpy as np
import math as mt

def computeAdjustedOperatorNBC(D2A, DD, tdex):
       # D2A is the operator to adjust
       # DD is the 1st derivative operator
       R = DD[tdex,tdex]
       dv = DD[tdex,:]; dv[tdex] = 0.0
       DA = (1.0 / R) * np.outer(DD[:,tdex], -dv)
              
       D2A[:,tdex] = 0.0
       D2A[tdex,:] = 0.0
       DOP = D2A + DA

       #DOP[tdex,:] = 0.0; 

       DOPC = numericalCleanUp(DOP)
       
       return DOPC

def computeAdjustedOperatorNBC_ends(D2A, DD):
       # D2A is the operator to adjust
       # DD is the 1st derivative operator
       
       R = (DD[0,0] * DD[-1,-1] - DD[0,-1] * DD[-1,0])
       
       lv = DD[0,:]; lv[0] = 0.0
       rv = DD[-1,:]; rv[-1] = 0.0
       
       V1 = (DD[-1,-1] / R) * (-lv) + (DD[0,-1] / R) * rv 
       DA1 = np.outer(DD[:,0], V1)
       
       V2 = (DD[-1,0] / R) * lv - (DD[0,0] / R) * (-rv) 
       DA2 = np.outer(DD[:,-1], V2)
       
       D2A[:,0] = 0.0; D2A[:,-1] = 0.0
       D2A[0,:] = 0.0; D2A[-1,:] = 0.0
       
       DOP = D2A + DA1 + DA2
       
       #DOP[0,:] = 0.0; DOP[:,0] = 0.0
       #DOP[-1,:] = 0.0; DOP[:,-1] = 0.0
       
       DOPC = numericalCleanUp(DOP)
       
       return DOPC

def computeAdjustedOperatorPeriodic(D2A):
       # D2A is the operator to adjust
       print('@@ PERIODIC BC @@')
       DOP = np.zeros(D2A.shape)
       NC = D2A.shape[1]
       
       # Copy over interior columns
       for jj in range(1,NC-1):
              DOP[:,jj] = D2A[:,jj] 
              
       # Couple the columns
       DOP[:,0] += D2A[:,-1]
       DOP[:,-1] += D2A[:,0]
       
       # Couple the rows and scale
       DOP[0,:] += D2A[-1,:]
       DOP[-1,:] += D2A[0,:]
       DOP[0,:] *= 0.5
       DOP[-1,:] *= 0.5
       
       return DOP

def numericalCleanUp(DDM):
       
       N = DDM.shape
       ZTOL = 1.0E-16
       DDMC = np.copy(DDM)
       # Clean up numerical zeros
       for ii in range(N[0]):
              for jj in range(N[1]):
                     if abs(DDM[ii,jj]) <= ZTOL:
                            DDMC[ii,jj] = 0.0
       return DDMC



# Computes Clamped Quintic Spline 1st derivative matrix
def computeQuinticSplineDerivativeMatrix(dom, DDM_BC):
       
       DM2 = DDM_BC.dot(DDM_BC)
       DM3 = DDM_BC.dot(DM2)
       DM4 = DDM_BC.dot(DM3)
       
       # Initialize matrix blocks
       N = len(dom)
       A = np.zeros((N,N)) # coefficients to 4th derivatives
       B = np.zeros((N,N)) # coefficients to RHS of 4th derivatives
       C = np.zeros((N,N)) # coefficients to 1st derivatives
       D = np.zeros((N,N)) # coefficients to additive part of 1st derivatives
       
       def computeIntegratalConstantMatrices(ii, x):
              
              x = np.array(x, dtype=np.longdouble)
              
              hp = abs(x[ii+1] - x[ii])
              hp1 = abs(x[ii+2] - x[ii+1])
              hm = abs(x[ii] - x[ii-1])
              
              V = np.zeros((9,9), dtype=np.longdouble)
              
              a0 = 0.5
              a1 = 1.0 / 6.0
              a2 = 1.0 / 24.0
              a3 = 1.0 / 120.0
              
              xim = x[ii-1]
              xi = x[ii]
              xi2 = a0*xi**2
              xip = x[ii+1]
              xip2 = a0*xip**2
              xir = x[ii+2]
              
              xqm = a1 * (xim**2 + xi * xim + xi**2)
              xqp = a1 * (xi**2 + xip * xi + xip**2)
              xqr = a1 * (xip**2 + xir * xip + xir**2)
              
              # A_j-1, B_j-1, C_j-1, A_j, B_j, C_j, A_j+1, B_j+1, C_j+1 
              V[0,:] = np.array([-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
              V[1,:] = np.array([-xi, -1.0, 0.0, xi, 1.0, 0.0, 0.0, 0.0, 0.0])
              V[2,:] = np.array([-xi2, -xi, -1.0, xi2, xi, 1.0, 0.0, 0.0, 0.0])
              V[3,:] = np.array([xqm, a0 * (xi + xim), 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
              V[4,:] = np.array([0.0, 0.0, 0.0, xqp, a0 * (xip + xi), 1.0, 0.0, 0.0, 0.0])
              V[5,:] = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, xqr, a0 * (xir + xip), 1.0])
              V[6,:] = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0])
              V[7,:] = np.array([0.0, 0.0, 0.0, xip, 1.0, 0.0, -xip, -1.0, 0.0])
              V[8,:] = np.array([0.0, 0.0, 0.0, xip2, xip, 1.0, -xip2, -xip, -1.0])
              
              rho = np.zeros((9,4), dtype=np.longdouble)
              rho[0,1] = a0 * (hp + hm)
              rho[1,1] = a1 * (hm**2 - hp**2)
              rho[2,1] = a2 * (hm**3 + hp**3)
              rho[3,0] = a3 * hm**3
              rho[3,1] = -a3 * hm**3
              rho[4,1] = a3 * hp**3
              rho[4,2] = -a3 * hp**3
              rho[5,2] = a3 * hp1**3
              rho[5,3] = -a3 * hp1**3
              rho[6,2] = -a0 * (hp + hp1)
              rho[7,2] = a1 * (hp1**2 - hp**2)
              rho[8,2] = -a2 * (hp**3 + hp1**3)
              
              eta = np.zeros((9,4), dtype=np.longdouble)
              eta[3,:] = 1.0 / hm * np.array([-1.0, 1.0, 0.0, 0.0])
              eta[4,:] = 1.0 / hp * np.array([0.0, -1.0, 1.0, 0.0])
              eta[5,:] = 1.0 / hp1 *np.array([0.0, 0.0, -1.0, 1.0])
              
              PLU = scl.lu_factor(V)
              VI = scl.lu_solve(PLU, np.eye(9))
              
              OM = VI.dot(rho)
              ET = VI.dot(eta)
              '''
              thisTol = 1.0E-15
              for ii in range(4):
                     OM[:,ii], es1 = sps.linalg.lgmres(sps.csc_matrix(V, dtype=np.float64), rho[:,ii], x0=OM[:,ii], M=VI, tol=thisTol, maxiter=10000)
                     ET[:,ii], es2 = sps.linalg.lgmres(sps.csc_matrix(V, dtype=np.float64), eta[:,ii], x0=ET[:,ii], M=VI, tol=thisTol, maxiter=10000)
              '''                     
              #OM = VIS.dot(rho)
              #ET = VIS.dot(eta)
                            
              return OM, ET
                     
       # Loop over each interior point in the irregular grid
       for ii in range(1,N-1):
              hp = abs(dom[ii+1] - dom[ii])
              hm = abs(dom[ii] - dom[ii-1])
              
              #%% LHS matrix
              a0 = 0.5
              a1 = 1.0 / 6.0
              a2 = 1.0 / 24.0
              hc = hp + hm              
              
              #'''
              if ii == 1:                     
                     # Compute and store left element integral coefficients
                     OM1 = np.zeros((3,4))
                     OM1[0,0] += a0 * hm
                     OM1[1,0] += -(a1 * hm**2 + a0 * hm * dom[0])
                     OM1[2,0] += (a2 * hm**3 + a1 * hm**2 * dom[0] + 0.25 * hm * dom[0]**2)
                     
                     ET1 = np.zeros((3,N))
                     ET1[0,:] = DM3[0,:]
                     ET1[1,:] = DM2[0,:] - dom[0] * DM3[0,:]
                     ET1[2,:] = DDM_BC[0,:] - dom[0] * DM2[0,:] + a0 * dom[0]**2 * DM3[0,:]
                     
                     # Compute the right EIC
                     OM2, ET2 = computeIntegratalConstantMatrices(ii, dom)
                     
                     # Assemble to the equation for Z
                     A[ii,ii] = -a0 * hc
                     A[ii,ii-1:ii+3] -= -(OM2[3,:])# * dom[ii] + OM2[1,:])
                     A[ii,0:4] -= +(OM1[0,:])# * dom[ii] + OM1[1,:])
                     
                     B[ii,ii-1:ii+3] += -(ET2[3,:])# * dom[ii] + ET2[1,:])
                     B[ii,:] += +ET1[0,:]
                     
                     # Compute the C matrix (coefficients to Z)
                     C[ii,ii] = -a2 * hp**3
                     C[ii,ii-1:ii+3] += a0 * dom[ii]**2 * OM2[3,:] + dom[ii] * OM2[4,:] + OM2[5,:]
                     
                     # Compute the D matrix (coefficients to Q)
                     D[ii,ii-1:ii+3] += a0 * dom[ii]**2 * ET2[3,:] + dom[ii] * ET2[4,:] + ET2[5,:]
                     
                     # Compute the C matrix (coefficients to Z)
                     C[0,0] = -a2 * hm**3
                     C[0,0] += a0 * dom[0]**2 * OM1[0,0] + dom[0] * OM1[1,0] + OM1[2,0]
                     
                     # Compute the D matrix (coefficients to Q)
                     D[0,:] += a0 * dom[0]**2 * ET1[0,:] + dom[0] * ET1[1,:] + ET1[2,:]
                     
              elif ii == N-2:
                     # Compute the left EIC
                     OM1, ET1 = computeIntegratalConstantMatrices(ii-1, dom)
                                          
                     # Compute and store right element integral coefficients
                     OM2 = np.zeros((3,4))
                     OM2[0,-1] += -a0 * hp
                     OM2[1,-1] += -(a1 * hp**2 - a0 * hp * dom[-1])
                     OM2[2,-1] += (-a2 * hp**3 + a1 * hp**2 * dom[-1] - 0.25 * hp * dom[-1]**2)
                     
                     ET2 = np.zeros((3,N))
                     ET2[0,:] = DM3[-1,:]
                     ET2[1,:] = DM2[-1,:] - dom[-1] * DM3[-1,:]
                     ET2[2,:] = DDM_BC[-1,:] - dom[-1] * DM2[-1,:] + a0 * dom[-1]**2 * DM3[-1,:]
                     
                     # Assemble to the equation for Z
                     A[ii,ii] = +a0 * hc
                     A[ii, N-4:N] -= -(OM2[0,:])# * dom[ii] + OM2[1,:])
                     A[ii, ii-2:ii+2] -= +(OM1[3,:])# * dom[ii] + OM1[1,:])
                     
                     B[ii,:] += -ET2[0,:]
                     B[ii, ii-2:ii+2] += +(ET1[3,:])# * dom[ii] + ET1[1,:])
                     
                     # Compute the C matrix (coefficients to Z)
                     C[ii,ii] = a2 * hm**3
                     C[ii,ii-2:ii+2] += a0 * dom[ii]**2 * OM1[3,:] + dom[ii] * OM1[4,:] + OM1[5,:]
                     
                     # Compute the D matrix (coefficients to Z)
                     D[ii,ii-2:ii+2] += a0 * dom[ii]**2 * ET1[3,:] + dom[ii] * ET1[4,:] + ET1[5,:]
                     
                     # Compute the C matrix (coefficients to Z)
                     C[-1,-1] = a2 * hp**3
                     C[-1,-1] += a0 * dom[-1]**2 * OM2[0,-1] + dom[-1] * OM2[1,-1] + OM2[2,-1]
                     
                     # Compute the D matrix (coefficients to Z)
                     D[-1,:] += a0 * dom[-1]**2 * ET2[0,:] + dom[-1] * ET2[1,:] + ET2[2,:]
              else:
                     # Compute adjacent EIC and assemble to internal equations for Z
                     OM1, ET1 = computeIntegratalConstantMatrices(ii-1, dom)
                     OM2, ET2 = computeIntegratalConstantMatrices(ii, dom)
                     
                     A[ii,ii] = -a0 * hc
                     A[ii, ii-1:ii+3] -= -(OM2[3,:])# * dom[ii] + OM2[1,:])
                     A[ii, ii-2:ii+2] -= +(OM1[3,:])# * dom[ii] + OM1[1,:])
                     
                     B[ii, ii-1:ii+3] += -(ET2[3,:])# * dom[ii] + ET2[1,:])
                     B[ii, ii-2:ii+2] += +(ET1[3,:])# * dom[ii] + ET1[1,:])
                     
                     # Compute the C matrix (coefficients to Z)
                     C[ii,ii] = -a2 * hp**3
                     C[ii,ii-1:ii+3] += a0 * dom[ii]**2 * OM2[3,:] + dom[ii] * OM2[4,:] + OM2[5,:]
                     
                     # Compute the D matrix (coefficients to Q)
                     D[ii,ii-1:ii+3] += a0 * dom[ii]**2 * ET2[3,:] + dom[ii] * ET2[4,:] + ET2[5,:]
              #'''
              
       # Left end natural S^(5) = 0
       #A[0,0] = 1.0 * dom[1] / (dom[1] - dom[0]) 
       #A[0,1] = -1.0 * dom[0] / (dom[1] - dom[0])
       #B[0,:] = np.zeros(N)
       # Left end known S^(4)
       A[0,0] = 1.0
       B[0,:] = DM4[0,:]
       
       # Right end natural S^(5) = 0
       #A[N-1,N-2] = 1.0 * dom[-1] / (dom[-1] - dom[-2])
       #A[N-1,N-1] = -1.0 * dom[-2] / (dom[-1] - dom[-2])
       #B[N-1,:] = np.zeros(N)
       # Right end known S^(4)
       A[-1,-1] = 1.0
       B[-1,:] = DM4[-1,:]
              
       # Compute the 4th derivative matrix
       PLU = scl.lu_factor(A)
       AIB = scl.lu_solve(PLU,B)
       
       # Compute the 1st derivative matrix
       DM1 = C.dot(AIB) + D
       
       #DM1[0,:] = DDM_BC[0,:]
       #DM1[-1,:] = DDM_BC[-1,:]
       DM1C = numericalCleanUp(DM1)
                     
       return DM1C, AIB

# Computes Cubic Spline 1st derivative matrix
def computeCubicSplineDerivativeMatrix(dom, isClamped, isEssential, \
                                       isLeftEssentialRightClamped, isLeftClampedRightEssential, DDM_BC):
       # Initialize matrix blocks
       N = len(dom)
       A = np.zeros((N,N)) # coefficients to 2nd derivatives
       B = np.zeros((N,N)) # coefficients to RHS of 2nd derivatives
       C = np.zeros((N,N)) # coefficients to 1st derivatives
       D = np.zeros((N,N)) # coefficients to additive part of 1st derivatives
       
       # Loop over each interior point in the irregular grid
       for ii in range(1,N-1):
              hp = abs(dom[ii+1] - dom[ii])
              hm = abs(dom[ii] - dom[ii-1])
              hc = abs(dom[ii+1] - dom[ii-1])
              
              A[ii,ii-1] = -1.0 / 6.0 * hm
              A[ii,ii] = -1.0 / 3.0 * hc
              A[ii,ii+1] = -1.0 / 6.0 * hp
              
              B[ii,ii-1] = -1.0 / hm
              B[ii,ii] = (1.0 / hm + 1.0 / hp)
              B[ii,ii+1] = -1.0 / hp
              
       for ii in range(1,N-1):
              hp = abs(dom[ii+1] - dom[ii])
              C[ii,ii] = -1.0 / 3.0 * hp
              C[ii,ii+1] = -1.0 / 6.0 * hp
              
              D[ii,ii] = -1.0 / hp
              D[ii,ii+1] = 1.0 / hp
              
       # Compute BC adjustments
       h0 = abs(dom[1] - dom[0])
       hn = abs(dom[N-1] - dom[N-2])
       
       if isClamped:
              # Left end
              A[0,0] = -2.0 * h0 / 3.0
              A[0,1] = h0 / 6.0
              
              # Use derivative by CFD to set boundary condition
              B[0,:] = DDM_BC[0,:]
              B[0,0] += 1.0 / h0
              B[0,1] -= 1.0 / h0
              
              # Right end
              A[N-1,N-2] = -hn / 6.0
              A[N-1,N-1] = 2.0 * hn / 3.0
              
              # Use derivative by CFD to set boundary condition
              B[N-1,:] = DDM_BC[N-1,:]
              B[N-1,N-2] += 1.0 / hn
              B[N-1,N-1] -= 1.0 / hn
              #'''
              
              # Compute the first derivative matrix
              AIB = np.linalg.solve(A, B)
              DDM = C.dot(AIB) + D
              
              # Adjust ends
              DDM[0,:] = DDM_BC[0,:]
              DDM[-1,:] = DDM_BC[-1,:]
              
       elif isEssential:
              # Left end
              A[0,0] = -1.0 / h0
              A[0,1] = 1.0 / h0
              B[0,:] = np.zeros(N)
              
              # Right end
              A[N-1,N-2] = 1.0 / hn
              A[N-1,N-1] = -1.0 / hn
              B[N-1,:] = np.zeros(N)
              
              # Compute the first derivative matrix
              AIB = np.linalg.solve(A, B)
              DDM = C.dot(AIB) + D
              
              # Adjust the ends
              DDM[0,:] = h0 / 6.0 * AIB[1,:] - 2.0 * h0 / 3.0 * AIB[0,:]
              DDM[0,0] -= 1.0 / h0
              DDM[0,1] += 1.0 / h0
              
              DDM[N-1,:] = -h0 / 6.0 * AIB[N-2,:] + 2.0 * hn / 3.0 * AIB[N-1,:]
              DDM[N-1,N-2] -= 1.0 / hn
              DDM[N-1,N-1] += 1.0 / hn
       
       elif isLeftEssentialRightClamped:
              # Left end
              A[0,0] = -1.0 / h0
              A[0,1] = 1.0 / h0
              B[0,:] = np.zeros(N)
              
              # Right end
              A[N-1,N-2] = -hn / 6.0
              A[N-1,N-1] = 2.0 * hn / 3.0
              
              # Use derivative by CFD to set boundary condition
              B[N-1,:] = DDM_BC[N-1,:]
              B[N-1,N-2] += 1.0 / hn
              B[N-1,N-1] -= 1.0 / hn
              
              # Compute the first derivative matrix
              AIB = np.linalg.solve(A, B)
              DDM = C.dot(AIB) + D
              #'''
              # Adjust ends
              DDM[0,:] = h0 / 6.0 * AIB[1,:] - 2.0 * h0 / 3.0 * AIB[0,:]
              DDM[0,0] -= 1.0 / h0
              DDM[0,1] += 1.0 / h0
              
              DDM[N-1,:] = 1.0 * DDM_BC[N-1,:]
              #'''
       elif isLeftClampedRightEssential:
              # Left end
              A[0,0] = -2.0 * h0 / 3.0
              A[0,1] = h0 / 6.0
              
              # Use derivative by CFD to set boundary condition
              B[0,:] = DDM_BC[0,:]
              B[0,0] += 1.0 / h0
              B[0,1] -= 1.0 / h0
              
              # Right end
              A[N-1,N-2] = 1.0 / hn
              A[N-1,N-1] = -1.0 / hn
              B[N-1,:] = np.zeros(N)
              
              # Compute the first derivative matrix
              AIB = np.linalg.solve(A, B)
              DDM = C.dot(AIB) + D
              #'''
              # Adjust ends
              DDM[0,:] = 1.0 * DDM_BC[0,:]
              
              DDM[N-1,:] = -h0 / 6.0 * AIB[N-2,:] + 2.0 * hn / 3.0 * AIB[N-1,:]
              DDM[N-1,N-2] -= 1.0 / hn
              DDM[N-1,N-1] += 1.0 / hn
              #'''
       else:
              # NATURAL cubic spline.
              AIB = np.zeros((N,N))
              AIB[1:N-1,1:N-1] = np.linalg.solve(A[1:N-1,1:N-1], B[1:N-1,1:N-1])
              DDM = C.dot(AIB) + D
              
              # Adjust the ends
              DDM[0,:] = h0 / 6.0 * AIB[1,:]
              DDM[0,0] -= 1.0 / h0
              DDM[0,1] += 1.0 / h0
              
              DDM[N-1,:] = -h0 / 6.0 * AIB[N-2,:]
              DDM[N-1,N-2] -= 1.0 / hn
              DDM[N-1,N-1] += 1.0 / hn      
              
       DDMC = numericalCleanUp(DDM)

       return DDMC, AIB

# Computes standard 4th order compact finite difference 1st derivative matrix
def computeCompactFiniteDiffDerivativeMatrix1(dom, order):
       
       end3 = False
       end6 = True
       
       # Initialize the left and right derivative matrices
       N = len(dom)
       LDM = np.zeros((N,N)) # tridiagonal
       RDM = np.zeros((N,N)) # centered difference
       
       def p2Matrix6(hm1, hp1, hp2, hp3, hp4, hp5, ND):
              CM = np.ones((ND,ND))
              c1 = hp1; c2 = hm1
              c3 = hp2 + hp1
              c5 = hp3 + hp2 + hp1
              c7 = hp4 + hp3 + hp2 + hp1
              #c9 = hp5 + hp4 + hp3 + hp2 + hp1

              # One sided boundary scheme (2 forward derivatives)
              CM[0,:] = [+c1, -c2, +c3, +c5, +c7, -2.0]
              CM[1,:] = [c1**2, +c2**2, +c3**2, +c5**2, +c7**2, -2.0 * (hp1 - hm1)]
              CM[2,:] = [c1**3, -c2**3, +c3**3, +c5**3, +c7**3, -3.0 * (hp1**2 + hm1**2)]
              CM[3,:] = [c1**4, +c2**4, +c3**4, +c5**4, +c7**4, -4.0 * (hp1**3 - hm1**3)]
              CM[4,:] = [c1**5, -c2**5, +c3**5, +c5**5, +c7**5, -5.0 * (hp1**4 + hm1**4)]
              CM[5,:] = [c1**6, +c2**6, +c3**6, +c5**6, +c7**6, -6.0 * (hp1**5 - hm1**5)]

              return CM
       
       def endMatrix6(hp1, hp2, hp3, hp4, hp5, ND):
              CM = np.ones((ND,ND))
              c1 = hp1
              c3 = hp2 + hp1
              c5 = hp3 + hp2 + hp1
              c7 = hp4 + hp3 + hp2 + hp1
              c9 = hp5 + hp4 + hp3 + hp2 + hp1

              # One sided boundary scheme (2 forward derivatives)
              CM[0,:] = [+c1, +c3, +c5, +c7, +c9, -1.0]
              CM[1,:] = [c1**2, +c3**2, +c5**2, +c7**2, +c9**2, -2.0 * c1]
              CM[2,:] = [c1**3, +c3**3, +c5**3, +c7**3, +c9**3, -3.0 * c1**2]
              CM[3,:] = [c1**4, +c3**4, +c5**4, +c7**4, +c9**4, -4.0 * c1**3]
              CM[4,:] = [c1**5, +c3**5, +c5**5, +c7**5, +c9**5, -5.0 * c1**4]
              CM[5,:] = [c1**6, +c3**6, +c5**6, +c7**6, +c9**6, -6.0 * c1**5]

              return CM
       
       def interiorMatrix10(hm3, hm2, hm1, hp1, hp2, hp3, ND):
              CM = np.ones((ND,ND))
              c1 = hp1; c2 = hm1
              c3 = hp2 + hp1; c4 = hm1 + hm2
              c5 = hp3 + hp2 + hp1; c6 = hm1 + hm2 + hm3

              #''' Pentadiagonal left, Septadiagonal right
              CM[0,:] = [c1, -c2, +c3, -c4, c5, -c6, -1.0, -1.0, -1.0, -1.0]
              CM[1,:] = [c1**2, +c2**2, +c3**2, +c4**2, +c5**2, +c6**2, +2.0 * c2, -2.0 * c1, +2.0 * c4, -2.0 * c3]
              CM[2,:] = [c1**3, -c2**3, +c3**3, -c4**3, +c5**3, -c6**3, -3.0 * c2**2, -3.0 * c1**2, -3.0 * c4**2, -3.0 * c3**2]
              CM[3,:] = [c1**4, +c2**4, +c3**4, +c4**4, +c5**4, +c6**4, +4.0 * c2**3, -4.0 * c1**3, +4.0 * c4**3, -4.0 * c3**3]
              CM[4,:] = [c1**5, -c2**5, +c3**5, -c4**5, +c5**5, -c6**5, -5.0 * c2**4, -5.0 * c1**4, -5.0 * c4**4, -5.0 * c3**4]
              CM[5,:] = [c1**6, +c2**6, +c3**6, +c4**6, +c5**6, +c6**6, +6.0 * c2**5, -6.0 * c1**5, +6.0 * c4**5, -6.0 * c3**5]
              CM[6,:] = [c1**7, -c2**7, +c3**7, -c4**7, +c5**7, -c6**7, -7.0 * c2**6, -7.0 * c1**6, -7.0 * c4**6, -7.0 * c3**6]
              CM[7,:] = [c1**8, +c2**8, +c3**8, +c4**8, +c5**8, +c6**8, +8.0 * c2**7, -8.0 * c1**7, +8.0 * c4**7, -8.0 * c3**7]
              CM[8,:] = [c1**9, -c2**9, +c3**9, -c4**9, +c5**9, -c6**9, -9.0 * c2**8, -9.0 * c1**8, -9.0 * c4**8, -9.0 * c3**8]
              CM[9,:] = [c1**10, +c2**10, +c3**10, +c4**10, +c5**10, +c6**10, +10.0 * c2**9, -10.0 * c1**9, +10.0 * c4**9, -10.0 * c3**9]
                            
              return CM
              
       for ii in range(1,N-1):
              # Get the metric weights
              hp1 = abs(dom[ii+1] - dom[ii])
              hm1 = abs(dom[ii] - dom[ii-1])
              
              if (order == 10 and ii in [2,N-3]) or \
                 (order == 6 and ii in range(2,N-2)):
                     hp2 = abs(dom[ii+2] - dom[ii+1])
                     hm2 = abs(dom[ii-1] - dom[ii-2])
              else:
                     hp2 = 0.0; hm2 = 0.0
              
              if (order == 10 and ii in range(3,N-3)):
                     hp2 = abs(dom[ii+2] - dom[ii+1])
                     hm2 = abs(dom[ii-1] - dom[ii-2])
                     hp3 = abs(dom[ii+3] - dom[ii+2])
                     hm3 = abs(dom[ii-3] - dom[ii-2])
              else:
                     hp3 = 0.0; hm3 = 0.0
              
              ND = 10
              CM10 = interiorMatrix10(hm3, hm2, hm1, hp1, hp2, hp3, ND)
              CMV = np.zeros(ND)
              CMV[0] = 1.0
       
              if (order == 4 and ii in range(1,N-1)):
                     
                     alpha = 0.25
                     beta = 0.25
                     
                     # Delete columns
                     ddex = [2, 3, 4, 5, 8, 9]
                     CM4 = np.delete(CM10, ddex, axis=1)
                                     
                     # Delete rows to 4th order
                     ddex = [4, 5, 6, 7, 8, 9]
                     CM4 = np.delete(CM4, ddex, axis=0)
                     CM4_V = np.delete(CMV, ddex, axis=0)
                     
                     # Constraint alpha = beta = 0.25
                     CM4[:,-2] += CM4[:,-1]
                     CM4_V -= alpha * CM4[:,-2]
                     CMS = CM4[0:-2,0:-2]
                     
                     PLU = scl.lu_factor(CMS)
                     CF = scl.lu_solve(PLU, CM4_V[0:-2])
                     
                     CFE = -np.sum(CF[0:2])
                     
                     # Write the right equation
                     RDM[ii,ii-1] = CF[1]
                     RDM[ii,ii] = CFE
                     RDM[ii,ii+1] = CF[0]
                     
                     # Write the left equation
                     LDM[ii,ii-1] = alpha
                     LDM[ii,ii] = 1.0
                     LDM[ii,ii+1] = beta
                     
              if (order > 4 and ii in [1,N-2]):
                     hm1 = dom[1] - dom[0]
                     hp1 = dom[2] - dom[1]
                     hp2 = dom[3] - dom[2]
                     hp3 = dom[4] - dom[3]
                     hp4 = dom[5] - dom[4]
                     hp5 = dom[6] - dom[5]
                     
                     CME = p2Matrix6(hm1, hp1, hp2, hp3, hp4, hp5, 6)
                                     
                     CME_V = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
                     
                     PLU = scl.lu_factor(CME)
                     CF = scl.lu_solve(PLU, CME_V)
                     
                     CFE = -np.sum(CF[0:-1])
                     
                     # Write the right equation
                     if ii == 1:
                            RDM[ii,ii-1] = CF[1]
                            RDM[ii,ii] = CFE
                            RDM[ii,ii+1] = CF[0]
                            RDM[ii,ii+2] = CF[2]
                            RDM[ii,ii+3] = CF[3]
                            RDM[ii,ii+4] = CF[4]
                     elif ii == N-2:
                            RDM[ii,ii+1] = -CF[1]
                            RDM[ii,ii] = -CFE
                            RDM[ii,ii-1] = -CF[0]
                            RDM[ii,ii-2] = -CF[2]
                            RDM[ii,ii-3] = -CF[3]
                            RDM[ii,ii-4] = -CF[4]
                     
                     # Write the left equation
                     LDM[ii,ii-1] = CF[-1]
                     LDM[ii,ii] = 1.0
                     LDM[ii,ii+1] = CF[-1]
              
              # Loop over each interior point in the irregular grid
              if (order == 10 and ii in [2,N-3]) or \
                 (order == 6 and ii in range(2,N-2)):   
                        
                     alpha = 1.0 / 3.0 
                     beta = 1.0 / 3.0
                     
                     # Delete columns
                     ddex = [4, 5, 8, 9]
                     CM6 = np.delete(CM10, ddex, axis=1) 
                     
                     # Delete rows to 4th order
                     ddex = [6, 7, 8, 9]
                     CM6 = np.delete(CM6, ddex, axis=0)
                     CM6_V = np.delete(CMV, ddex, axis=0)
                     
                     # Constraint alpha = beta = 1/3
                     CM6[:,-2] += CM6[:,-1]
                     CM6_V -= alpha * CM6[:,-2]
                     CMS = CM6[0:-2,0:-2]
                     
                     PLU = scl.lu_factor(CMS)
                     CF = scl.lu_solve(PLU, CM6_V[0:-2])
                     
                     CFE = -np.sum(CF[0:4])
                     
                     # Write the right equation
                     RDM[ii,ii-2] = CF[3]
                     RDM[ii,ii-1] = CF[1]
                     RDM[ii,ii] = CFE
                     RDM[ii,ii+1] = CF[0]
                     RDM[ii,ii+2] = CF[2]
                     
                     # Write the left equation
                     LDM[ii,ii-1] = alpha
                     LDM[ii,ii] = 1.0
                     LDM[ii,ii+1] = beta
                     
              # Loop over each interior point in the irregular grid
              if (order == 10 and ii in range(3,N-3)):
                     
                     alpha = 0.5
                     beta = 0.5
                     theta = 0.05
                     rho = 0.05
                     
                     CMI = np.copy(CM10)
                     CMI[:,-2] += CMI[:,-1]
                     CMI[:,-4] += CMI[:,-3]
                     CM10_V = CMV - (alpha * CMI[:,-4] + theta * CMI[:,-2])
                     
                     sdex = np.array([0, 1, 2, 3, 4, 5])
                     CMS = CMI[np.ix_(sdex,sdex)]
                     
                     PLU = scl.lu_factor(CMS)
                     CF = scl.lu_solve(PLU, CM10_V[sdex])
                     CFE = -np.sum(CF[0:6])
                     
                     # Write the right equation
                     RDM[ii,ii-3] = CF[5]
                     RDM[ii,ii-2] = CF[3]
                     RDM[ii,ii-1] = CF[1]
                     RDM[ii,ii] = CFE
                     RDM[ii,ii+1] = CF[0]
                     RDM[ii,ii+2] = CF[2]
                     RDM[ii,ii+3] = CF[4]
                     
                     # Write the left equation
                     LDM[ii,ii-2] = theta
                     LDM[ii,ii-1] = alpha
                     LDM[ii,ii] = 1.0
                     LDM[ii,ii+1] = beta
                     LDM[ii,ii+2] = rho
       
       if end6: 
              # Coefficients for 6th order compact one-sided schemes
              hp1 = dom[1] - dom[0]
              hp2 = dom[2] - dom[1]
              hp3 = dom[3] - dom[2]
              hp4 = dom[4] - dom[3]
              hp5 = dom[5] - dom[4]
       
              # Compute the stencil coefficients
              ND = 6
              CME = endMatrix6(hp1, hp2, hp3, hp4, hp5, 6)
              CME_V = np.zeros(ND)
              CME_V[0] = 1.0
              
              PLU = scl.lu_factor(CME)
              CF_F = scl.lu_solve(PLU, CME_V)
              beta = CF_F[-1]
              
              LDM[0,0] = 1.0
              LDM[0,1] = beta
              RDM[0,0] = -np.sum(CF_F[0:-1])
              RDM[0,1] = CF_F[0]
              RDM[0,2] = CF_F[1]
              RDM[0,3] = CF_F[2]
              RDM[0,4] = CF_F[3]
              RDM[0,5] = CF_F[4]
              
              LDM[N-1,N-1] = 1.0
              LDM[N-1,N-2] = beta
              RDM[N-1,N-1] = np.sum(CF_F[0:-1])
              RDM[N-1,N-2] = -CF_F[0]
              RDM[N-1,N-3] = -CF_F[1]
              RDM[N-1,N-4] = -CF_F[2]
              RDM[N-1,N-5] = -CF_F[3]
              RDM[N-1,N-6] = -CF_F[4]
       
       if end3:
              hp2 = dom[1] - dom[0]
              hp3 = dom[2] - dom[1]
              LDM[0,0] = 1.0
              LDM[0,1] = (hp2 + hp3) / hp3
              RDM[0,0] = -(3.0 * hp2 + 2.0 * hp3) / (hp2 * (hp2 + hp3))
              RDM[0,1] = (hp2 + hp3) * (2.0 * hp3 - hp2) / (hp2 * hp3**2)
              RDM[0,2] = (hp2**2) / (hp3**2 * (hp2 + hp3))
              '''
              '''
              hp2 = dom[N-2] - dom[N-1]
              hp3 = dom[N-3] - dom[N-2]
              LDM[N-1,N-1] = 1.0
              LDM[N-1,N-2] = (hp2 + hp3) / hp3
              RDM[N-1,N-1] = -(3.0 * hp2 + 2.0 * hp3) / (hp2 * (hp2 + hp3))
              RDM[N-1,N-2] = (hp2 + hp3) * (2.0 * hp3 - hp2) / (hp2 * hp3**2)
              RDM[N-1,N-3] = (hp2**2) / (hp3**2 * (hp2 + hp3))
              
       # Get the derivative matrix
       DDM = np.linalg.solve(LDM, RDM)
       
       DDMA = numericalCleanUp(DDM)
       
       return DDMA
       
       # Get data from DIMS
       ZH = DIMS[2]
       NZ = DIMS[4]
       
       # Initialize grid and make column vector
       xi, wcp = cheblb(NZ)
       
       b = 2.0 / ZH
   
       # Get the Chebyshev transformation matrix
       CT = chebpolym(NZ+1, -xi)
   
       # Make a diagonal matrix of weights
       W = np.diag(wcp)
   
       # Compute scaling for the forward transform
       S = np.eye(NZ+1)
   
       for ii in range(NZ):
              temp = W.dot(CT[:,ii])
              temp = ((CT[:,ii]).T).dot(temp)
              S[ii,ii] = temp ** (-1)

       S[NZ,NZ] = 1.0 / mt.pi
   
       # Compute the spectral derivative coefficients
       SDIFF = np.zeros((NZ+1,NZ+1))
       SDIFF[NZ-1,NZ] = 2.0 * NZ
   
       for ii in reversed(range(NZ - 1)):
              A = 2.0 * (ii + 1)
              B = 1.0
              if ii > 0:
                     c = 1.0
              else:
                     c = 2.0
            
              SDIFF[ii,:] = B / c * SDIFF[ii+2,:]
              SDIFF[ii,ii+1] = A / c
    
       # Chebyshev spectral transform in matrix form
       temp = CT.dot(W)
       STR_C = S.dot(temp)
       # Chebyshev spatial derivative based on spectral differentiation
       # Domain scale factor included here
       temp = (CT.T).dot(SDIFF)
       DDM = -b * temp.dot(STR_C)
       
       DDMC = numericalCleanUp(DDM)
              
       #print(xi)
       #print(DDM[0,0], -(2.0 * NZ**2 + 1) / 3.0 / ZH)
       #print(DDM[-1,-1], (2.0 * NZ**2 + 1) / 3.0 / ZH)

       return DDMC, STR_C

def computeFourierDerivativeMatrix(DIMS):
       
       # Get data from DIMS
       L1 = DIMS[0]
       L2 = DIMS[1]
       NX = DIMS[3]
       
       kxf = (2*mt.pi/abs(L2 - L1)) * np.fft.fftfreq(NX+1) * (NX+1)
       KDM = np.diag(kxf, k=0)
       DFT = np.fft.fft(np.eye(NX+1), axis=0)
       DDM = np.fft.ifft(1j * KDM.dot(DFT), axis=0)
       
       DDMC = numericalCleanUp(DDM)
       
       return DDMC, DFT