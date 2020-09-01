import matplotlib.pyplot as plt
import numpy as np
import csv

kB = 1.38e-23

def InitialConditions(x, N):
  ##--Initial Conditions----
  xmax = max(x)
  xmin = min(x)
  rho = []
  vel = []
  Pressure = []
  engyDens = []

  for i in range(1, N+1):
    rho.append( np.sin(2*np.pi*x[i] / (xmax - xmin)))#rho.append(np.exp(-beta * (i - 0.5)**2))
    engyDens.append(np.exp(-(x[N/2] - x[i])**2)) #Energy Density Gaussian
    if (i < 60 ):
#      rho.append(0.)
      vel.append(-0.5)
    elif (i >= 60):# < 85):
#      rho.append(0.5)
      vel.append(0.5) 
  #  else:
  #    rho.append(0.)
#    vel.append(1.)  

  rho_original = rho
  Pressure = vel #initial Pressure

  return rho, vel, engyDens, Pressure, rho_original

def Temp(Ndens, engyDens, rho, vel, N):
  Temp = []
  for i in range(1,N+1):
    Temp.append((2./5.) * (engyDens[i] - 0.5 * rho[i] * vel[i] * vel[i]) / (Ndens*kB))
  return Temp

def PowerDeposition(N,dt):
  Q = []
  for i in range(1,N+1):
    Q.append(i*dt*100)
  return Q

def Upwind(a, dx, N):
  b = []
  for i in range(1,N+1):
    b.append((a[i]-a[i-1])/dx)
  return b


def fill_ghosts(q, ng, N):
  for i in range(ng-1):
    #left Boundary
    q.insert(0, q[N-1])
#    q[0] = q[N]
    #right Boundary 
    q.insert(N+1, q[0])
#    q[N-1+i] = q[ng+i]
  return q


def states(a, dx, dt, vel, slope_type, N):
  if(slope_type == 'SuperBee'):
    C = 0.
    D = 0.
    Ql = []
    Qr = []
    for i in range(1,N+1):
      A = minmod1(C, D, a, dx, i)
      B = minmod2(C, D, a, dx, i)
      slope = MaxMod(A, B, i)
#      print(slope)     
 
      Ql.append(a[i-2] + 0.5 * (dx - vel[i]*dt) * slope)
      Qr.append(a[i-1] - 0.5 * (dx - vel[i]*dt) * slope)

#    for i in range(N):
  return Qr, Ql


def minmod1(C, D, a, dx, i):
  C = (a[i+1] - a[i])/dx
  D = 2.*(a[i]-a[i-1])/dx
  if (abs(C) < abs(D) and C*D > 0.):
    return C
  elif (abs(D) < abs(C) and C*D > 0.):
    return D
  else:
    return 0.


def minmod2(C, D, a, dx, i):
  C = (a[i] - a[i-1])/dx
  D = 2.*(a[i+1]-a[i])/dx

  if(abs(C) < abs(D) and C*D > 0.):
    return C
  elif(abs(D) < abs(C) and C*D > 0.):
    return D
  else:
    return 0.

def MaxMod(A, B, i):
  if(abs(A) < abs(B) and A*B > 0.):
    return A
  elif(abs(B) < abs(A) and A*B > 0.):
    return B
  else:
    return 0.

def Riemann(Ql, Qr, vel, N):
  flux = []
  for i in range(N):
    if (vel[i] > 0.):
      flux.append(Ql[i]*vel[i+1])
    else:
      flux.append(Qr[i]*vel[i])

  return flux
