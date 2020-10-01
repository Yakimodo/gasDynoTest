import matplotlib.pyplot as plt
import numpy as np
import csv

kB = 1.38e-23

def initialConditions(sim_type, sim_details, xc, N):
  ##--Initial Conditions----
  xmax = max(xc)
  xmin = min(xc)
  rho = []
  u = []
  Pressure = []
  engyDens = []
  q = np.empty((N+2,2))
  s = []
  K = []
  c_sound = []
  c_water = []
  Z = []

  for i in range(N):

    if (i <= 0.5*N):

      if (sim_type == 'ShallowWater' and sim_details == 'SmallDisturbance'):
        Pressure.append(1+np.exp(-np.abs((xc[i]-max(xc)//2.))**2))#h/P #(3.) #np.sqrt(abs(1.-(xc[i]-3)**2))*(xc[i]<4.)*(xc[i]>2.)) ##Interface--
        u.append(0.)
      elif (sim_type == 'ShallowWater' and sim_details == 'DamBreak'):
        Pressure.append(3.)
        u.append(0.)

      rho.append(1.) #(1-homoInter)__(1-heterInter)___(3---shockTube)
      K.append(1.)
      c_sound.append(np.sqrt(K[i]/rho[i]))#Sound
      Z.append(rho[i]*c_sound[i])

      if (sim_type == 'Sound'):
        Pressure.append(np.sqrt(abs(1.-(xc[i]-3)**2))*(xc[i]<4.)*(xc[i]>2.))
        u.append(Pressure[i]/Z[i])

#      Pressure.append(3.)

      c_water.append(np.sqrt(Pressure[i])) #ShallowWater

    elif (i > 0.5*N):

      if (sim_type == 'ShallowWater' and sim_details == 'SmallDisturbance'):
        Pressure.append(1+np.exp(-np.abs((xc[i]-max(xc)//2.))**2))#h/P #(3.) #np.sqrt(abs(1.-(xc[i]-3)**2))*(xc[i]<4.)*(xc[i]>2.)) ##Interface--
        u.append(0.)
      elif (sim_type == 'ShallowWater' and sim_details == 'DamBreak'):
        Pressure.append(1.)
        u.append(0.)

      rho.append(4.) #(2-homoInter)__(4-heterInter)___(1---shockTube)
      K.append(1.)
      c_sound.append(np.sqrt(K[i]/rho[i])) #Sound
      Z.append(rho[i]*c_sound[i]) #Zl = Zr for constant coefficient case

      if (sim_type == 'Sound'):
        Pressure.append(np.sqrt(abs(1.-(xc[i]-3)**2))*(xc[i]<4.)*(xc[i]>2.))
        u.append(Pressure[i]/Z[i])
#      Pressure.append(1.)

      c_water.append(np.sqrt(Pressure[i])) #ShallowWater

  for j in range(N):
      q[j+1,0] = Pressure[j]
      if (sim_type == 'Sound'):
        q[j+1,1] = u[j]
      if (sim_type == 'ShallowWater'):
        q[j+1,1] = u[j]*Pressure[j]

  rho_original = rho
#    rho.append( np.sin(2*np.pi*xc[i] / (xmax - xcmin)))#rho.append(np.exp(-beta * (i - 0.5)**2))
#    rho.append(np.sin(xc[i]))
#    rho.append(np.exp(-(xc[N/2] - xc[i])**2)) #Energy Density Gaussian
#    engyDens.append(np.exp(-(xc[N/2] - xc[i])**2)) #Energy Density Gaussian

  return u, Pressure, q, Z, K, rho

def Mat_ghostCells(q, N, BC):
  if(BC == 'extrap'):
    q[0,0] = q[1,0]
    q[N+1,0] = q[N,0]
    q[0,1] = q[1,1]
    q[N+1,1] = q[N,1]
  return q

def fill_ghosts(q, N, Type):
#  for i in range(ng-1):
  if (Type == 'Periodic'):
    #left Boundary
    q.insert(0, q[-1])
    #right Boundary 
    q.insert(len(q), q[1])
  
  if (Type == 'extrap'):
    q.insert(0, q[0])
    q.insert(len(q), q[-1])
  return q


def JumpSplit(q, N, i): 

  dq1 = 0.
  dq2 = 0.
#  for i in range(1,N+1):
  dq1 = q[i,0] - q[i-1,0]
  dq2 = q[i,1] - q[i-1,1]   

  return dq1,dq2




def Temp(Ndens, engyDens, rho, u, N):
  Temp = []
  for i in range(1,N+1):
    Temp.append((2./5.) * (engyDens[i] - 0.5 * rho[i] * u[i] * u[i]) / (Ndens*kB))
  return Temp

def PowerDeposition(N,t):
  Q = []
  for i in range(1,N+1):
    Q.append(i*t*10)
  return Q

def Upwind(a, dxc, N, BC):
  b = []
  for i in range(1,N+1):
    b.append((a[i]-a[i-1])/dxc)
  b = fill_ghosts(b, 2, N, BC)
  return b

def CentDiff(a, dxc, N, BC):
  b = []
  for i in range(1,N+1):
    b.append((a[i+1]-a[i-1])/(2*dxc))
  b = fill_ghosts(b, 2, N, BC)
  return b

def states(a, dxc, dt, u, slope_type, N, BC):
  C = 0.
  D = 0.
  E = 0.
  Ql = []
  Qr = []
  slope = []
  if(slope_type == 'SuperBee'):
    for i in range(1, N+1):
      A = minmod1(C, D, a, dxc, i)
      B = minmod2(C, D, a, dxc, i)
      slope.append(MaxMod(A, B))

  if(slope_type == 'MC Limiter'):
    for i in range(1,N+1):
      slope.append(minmodMC(C,D,E, a, dxc, i)) 
  slope = fill_ghosts(slope, 2, N, BC)
  for i in range(1, N+1): 
    Ql.append(a[i-1] + 0.5 * (dxc - u[i-1]*dt) * slope[i-1])
    Qr.append(a[i] - 0.5 * u[i-1] * (dxc - u[i-1]*dt) * slope[i])
  print(len(Ql))
  return Qr, Ql


def minmod1(C, D, a, dxc, i):
  C = (a[i+1] - a[i])/dxc
  D = 2.*(a[i]-a[i-1])/dxc
  
  if (abs(C) < abs(D) and C*D > 0.):
    return C
  elif (abs(D) < abs(C) and C*D > 0.):
    return D
  else:
    return 0.


def minmod2(C, D, a, dxc, i):
  C = (a[i] - a[i-1])/dxc
  D = 2.*(a[i+1]-a[i])/dxc

  if(abs(C) < abs(D) and C*D > 0.):
    return C
  elif(abs(D) < abs(C) and C*D > 0.):
    return D
  else:
    return 0.

def minmodMC(C, D, E, a, dxc, i):
  C = (a[i+1] - a[i-1])/(2.*dxc) 
  D = 2.*(a[i]-a[i-1])/dxc
  E = 2.*(a[i+1]-a[i])/dxc

  if(abs(C) < abs(D) and abs(C) < abs(E)  and C*D*E > 0.):
    return C
  elif(abs(D) < abs(C) and abs(D) < abs(E) and C*D*E > 0.):
    return D
  elif(abs(E) < abs(C) and abs(E) < abs(E) and C*D*E > 0.):
    return E
  else:
    return 0.


def MaxMod(A, B):
  if(abs(A) > abs(B) and A*B > 0.):
    return A
  elif(abs(B) > abs(A) and A*B > 0.):
    return B
  else:
    return 0.

def Riemann(Ql, Qr, u, N):
  flux = []
  for i in range(N):
    if (u[i] > 0.):
      flux.append(Ql[i]*u[i])
    else:
      flux.append(Qr[i]*u[i+1])

  return flux

def VelChange(u, N, t):
  a = []
  for i in range(N):
    a.append(u[i]*1.1)
  return a

def Porta_plotty2(X, y, y_original, MOVIE, Title, pauseRate, figNum, t, dt, NumPlots):
  if (MOVIE == 1):
    plt.figure(figNum)
    plt.clf()
    plt.title(Title+ ' ' +'Density Time $t$ = ' + str(t*dt) + ' ns')
    plt.plot(X[:], y[:]) #, label = str((n-1) * 0.1) + ' ns')
#  plt.plot(x[1:-1], rho_new1) #, label = str((n-1) * 0.1) + ' ns'
    plt.plot(X, y_original[1:-1])
#   plt.axhline(y=1, color='r', linestyle='--')
#  plt.axis((0, 2*np.pi, -1, 1)) 
    plt.pause(pauseRate)
    plt.legend()

  if (MOVIE == 0):
    plt.figure(figNum)
    plt.title(Title + ' ' +'Density Time $t$ = ' + str(t*dt) + ' ns')
    plt.plot(X[:], y[:]) #, label = str((n-1) * 0.1) + ' ns')
#  plt.plot(x[1:-1], rho_new1) #, label = str((n-1) * 0.1) + ' ns'
#  plt.plot(x, rho_original[:])
#   plt.axhline(y=1, color='r', linestyle='--')
#  plt.axis((0, 2*np.pi, -1, 1)) 
#    plt.pause(pauseRate)
    plt.legend()

