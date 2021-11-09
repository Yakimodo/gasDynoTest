#======function_list.py contains:
#---initialConditions(sim_type, sim_details, xo, xc, N, q):
#---Mat_ghostCells(q, N, BC)
#---fill_ghosts(q, N, Type)
#---JumpSplit(q, N, sim_type):
#---Temp(Ndens, engyDens, rho, u, N)
#---PowerDeposition(N,t)
#---Upwind(a, dxc, N, BC)
#---splitValues(qMat)
#---CentDiff(a, dxc, N, BC)
#---states(a, dxc, dt, u, slope_type, N, BC):
#---minmod1(C, D, a, dxc, i)
#---minmod2(C, D, a, dxc, i)
#---minmodMC(C, D, E, a, dxc, i)
#---MaxMod(A, B)
#---Riemann(Ql, Qr, u, N)
#---VelChange(u, N, t)
#---Porta_plotty2(X, y, y_original, MOVIE, Title, pauseRate, figNum, t, dt, NumPlots):	
#---
#---
#---
#---



import numpy as np
import matplotlib.pyplot as plt

kB = 1.38e-23
gamma = 1.4 #air == 7/5
beta = 0.2

def initialConditions(sim_type, sim_details, xc, xo, N, q):
  ##--Initial Conditions----
  xmax = max(xc)
  xmin = min(xc)
  rho = []
  u = []
  Pressure = []
  height = []
  engyDens = []
#  q = np.empty((3, N+2))
  s = []
  K = []
  c_sound = []
  c_water = []
  c_gas = []
  Z = []
  gamma1 = gamma - 1  
  g = 1. # gravity for shallow water eqs.


######________Shallow Water___________
  if (sim_type == 'ShallowWater'):
    for i in range(N):

      if (sim_details == 'SmallDisturbance'):
        height.append(1.+0.4*np.exp(-np.abs((xc[i]-xc[N//2]))**2 / beta))
        u.append(0.)

      elif (sim_details == 'DamBreak'):
        if (i <= N/2):
          height.append(1.5)
          u.append(0.)
        elif (i > N/2):
          height.append(1.)
          u.append(0.)

      q[0, i+1] = height[i]
      q[1, i+1] = u[i]*height[i]
    if q.shape[0] == 3:
      q[2, :] = 0.


######_________Sound_________________
  if (sim_type == 'Sound'):
    for i in range(N):
      if (i <= 0.5*(N)):
        Pressure.append(np.sqrt(abs(1.-(xc[i]-3)**2))*(xc[i]<4.)*(xc[i]>2.))
        u.append(Pressure[i]/Z[i])
        rho.append(1.) #(1-homoInter)__(1-heterInter)___(3---shockTube)
        K.append(1.)
        c_sound.append(np.sqrt(K[i]/rho[i])) #Sound
        Z.append(rho[i]*c_sound[i])

      elif (i > 0.5*(N)):
        Pressure.append(np.sqrt(abs(1.-(xc[i]-3)**2))*(xc[i]<4.)*(xc[i]>2.))
        u.append(Pressure[i]/Z[i])
        rho.append(4.) #(2-homoInter)__(4-heterInter)___(1---shockTube)
        K.append(1.)
        c_sound.append(np.sqrt(K[i]/rho[i])) 
        Z.append(rho[i]*c_sound[i]) #Zl = Zr for constant coefficient case
#      Pressure.append(1.)

      q[0, j+1] = Pressure[j]
      q[1, j+1] = u[j]
      q[2, j+1] = Z[j]

######____________________________________############
######_________GAS DYNAMICS_______________############
######____________________________________############

  if (sim_type == 'gasDyno'):

    if (sim_details == 'shockTube'):
      print('====SHOCK TUBE=====')
      for j in range(N):
#        if (j <= (N-1)/2.): #LEFT STATE
        if (xc[j] <= xo*max(xc)): #TORO LEFT STATE 
          rho.append(1.)
          u.append(.75) #0 or 0.75
          Pressure.append(1.) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] *  u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))
  
#        elif (j > (N-1)/2.): #RIGHT STATE
        elif (xc[j] > xo*max(xc)): #TORO RIGHT STATE 
          rho.append(0.125)
          u.append(0.)
          Pressure.append(0.1) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] * u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))

        q[0, j+1] = rho[j]
        q[1, j+1] = rho[j]*u[j] #u[j]
        q[2, j+1] = engyDens[j] #Pressure[j]

    if (sim_details == 'ToroPg151_Test2'):
      print('====ToroPg151_Test2=====')
      for j in range(N):
        if (j <= (N-1)/2.): #LEFT STATE
          rho.append(1.)
          u.append(-2.)
          Pressure.append(0.4) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] *  u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))
  
        elif (j > (N-1)/2.): #RIGHT STATE
          rho.append(1.)
          u.append(2.)
          Pressure.append(0.4) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] * u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))

        q[0, j+1] = rho[j]
        q[1, j+1] = rho[j]*u[j] #u[j]
        q[2, j+1] = engyDens[j] #Pressure[j]

    if (sim_details == 'ToroPg151_Test3'):
      print('====ToroPg151_Test3=====')
      for j in range(N):
        if (xc[j] <= xo*max(xc)): #LEFT STATE
          rho.append(1.)
          u.append(0.)
          Pressure.append(1000.) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] *  u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))
  
        elif (xc[j] > xo*max(xc)): #RIGHT STATE
          rho.append(1.)
          u.append(0.)
          Pressure.append(0.01) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] * u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))

        q[0, j+1] = rho[j]
        q[1, j+1] = rho[j]*u[j] #u[j]
        q[2, j+1] = engyDens[j] #Pressure[j]

    if (sim_details == 'ToroPg151_Test4'):
      print('====ToroPg151_Test4=====')
      for j in range(N):
        if (xc[j] <= xo*max(xc)): #LEFT STATE
          rho.append(5.9924)
          u.append(19.5975)
          Pressure.append(460.894) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] *  u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))
  
        elif (xc[j] > xo*max(xc)): #RIGHT STATE
          rho.append(5.99242)
          u.append(-6.19633)
          Pressure.append(46.0950) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] * u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))

        q[0, j+1] = rho[j]
        q[1, j+1] = rho[j]*u[j] #u[j]
        q[2, j+1] = engyDens[j] #Pressure[j]

    if (sim_details == 'ToroPg151_Test5'):
      print('====ToroPg151_Test5=====')
      for j in range(N):
        if (j <= (N-1)/2.): #LEFT STATE
          rho.append(5.99924)
          u.append(19.5975)
          Pressure.append(460.894) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] *  u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))
  
        elif (j > (N-1)/2.): #RIGHT STATE
          rho.append(5.99242)
          u.append(-6.19633)
          Pressure.append(46.0950) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] * u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))

        q[0, j+1] = rho[j]
        q[1, j+1] = rho[j]*u[j] #u[j]
        q[2, j+1] = engyDens[j] #Pressure[j]

    if (sim_details == 'BalbasICs'):
      print('====Balbas=====')
      for j in range(N):
        if (j <= (N-1)/2.): #LEFT STATE
          u.append(1.206)
          rho.append(0.5323)
          Pressure.append(0.3) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] *  u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))
  
        elif (j > (N-1)/2.): #RIGHT STATE
          u.append(1.206)
          rho.append(0.138)
          Pressure.append(0.029) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] * u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))

        q[0, j+1] = rho[j]
        q[1, j+1] = rho[j]*u[j] #u[j]
        q[2, j+1] = engyDens[j] #Pressure[j]

    if (sim_details == 'twoShock'):
      print('====TWO SHOCK====')
          
      for j in range(N):
        if (j <= (N-1)/2.): #LEFT STATE
          u.append(3.)
          rho.append(1.)
          Pressure.append(1.) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] *  u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))
  
        elif (j > (N-1)/2.): #RIGHT STATE
          u.append(1.)
          rho.append(2.)
          Pressure.append(1.) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] * u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))

        q[0, j+1] = rho[j]
        q[1, j+1] = rho[j]*u[j] #u[j]
        q[2, j+1] = engyDens[j] #Pressure[j]

    if (sim_details == 'flat'):
      print('====flat=====')
      for j in range(N):
        if (j <= (N-1)/2.): #LEFT STATE
          u.append(0.)
          rho.append(1.)
          Pressure.append(0.) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] *  u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))
  
        elif (j > (N-1)/2.): #RIGHT STATE
          u.append(0.)
          rho.append(1.)
          Pressure.append(0.) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] * u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))

        q[0, j+1] = rho[j]
        q[1, j+1] = rho[j]*u[j] #u[j]
        q[2, j+1] = engyDens[j] #Pressure[j]

  return q, u, Pressure#,, Z, K, rho

def Mat_ghostCells(q, N, BC):
  if(BC == 'extrap'):
    q[0,0] = q[0,1]
    q[0,-1] = q[0,-2]
    q[1,0] = q[1,1]
    q[1,-1] = q[1,-2]
    if q.shape[0] == 3:
      q[2,0] = q[2,1]
      q[2,-1] = q[2,-2]
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


def JumpSplit(q, N, sim_type): 
  dq1 = 0.
  dq2 = 0.
  dq3 = 0.
#  for i in range(1,N+1):
  dq1 = q[0,1:] - q[0,:-1] #q[0,i] - q[0,i-1]
  dq2 = q[1,1:] - q[1,:-1] #q[1,i] - q[1,i-1]
  if (sim_type == 'gasDyno'):   
    dq3 = q[2,1:] - q[2,:-1] #q[2,i] - q[2,i-1]   
  else:
    dq3 = 0.

  return dq1,dq2,dq3




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

def splitValues(qMat):
  ql = np.empty((qMat.shape[0], qMat.shape[1]-1))
  qr = np.empty((qMat.shape[0], qMat.shape[1]-1))
  for i in range(qMat.shape[0]):
    for j in range(qMat.shape[1]):
      if (j < qMat.shape[1]-1):
        ql[i,j] = qMat[i,j] - qMat[i,j-1]
      if (j > 0):
        qr[i,j] = qMat[i,j+1] - qMat[i,j]
  return ql, qr

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


#def FEMethod():
#  un+1 = un + dt*(ui+1 - ui)



#for RHS of charged species continuity equation dn/dt = dj/dx + kinetic Source Terms  ==>
#                                  DISCRETIZED:  n(t+1) = n(t) + dt/dx * ( j(k+1/2) - j(k-1/2) ) + dt * (Kinetic Source terms)
def SG(x, D, mobi, n, E):

  for k in range(len(x)):
    hk = 0.5*(x[k] + x[k+1]) 
    Dk = 0.5*(D[k] + D[k+1])
    mu = 0.5*(mobi[k] + mobi[k+1])
    Ek = 0.5*(E[k]+E[k+1])
    v_dr = mu*Ek
    
    alpha = v_dr*hk / Dk
  
    Io = (np.exp(alpha) - 1 ) / alpha

    flux_density.append( Dk / hk * (n[k] - np.exp(alpha)*n[k+1]) / Io )

  return flux_density




def SOR(s, z, N, CYL, phi):

  ds = s[i+1] - s[i]
  dz = z[i+1] - z[i]
  RT = 1e-6
  AT = 0.1

  if (CYL == 0):
    for i in range(N):
      r = 0
      for m in range(MaxIts):
        while r == 0:
          w = 2. - 2.*np.pi*ds
          phi_new[i] = w/2. * (phi_new[i-1] + phi_old[i+1] + ds**2*delta) + (1-w)*phi_old
          W[i] = 1. / (RT*np.abs(phi_new[i] + AT))
          WRMS = np.sqrt( 1./N * np.sum((W[i] * np.abs(phi_new[i] - phi_old[i]))**2))   
          if (WRMS < 1.):
            r = 1          
          

    
#  if (CYL == 1):
#    if(s != 0 or s != N or z != 0 or z != N):
#      EE = - ( 2/ds**2 + 2/dz**2)
#      A = (1/dz**2)/EE
#      B = A
#      C = (1

  return phi_new


##=== E = grad(phi)
#def calcField(phi_new):
#  E =  
  




