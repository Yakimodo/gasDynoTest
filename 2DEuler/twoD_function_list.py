import numpy as np
import matplotlib.pyplot as plt

kB = 1.38e-23
gamma = 1.4 #air == 7/5
beta = 0.2

def initialConditions(sim_type, sim_details, rc, zc, M, N, q):
  ##--Initial Conditions----
  rmax = max(rc)
  rmin = min(rc)
  zmax = max(zc)
  zmin = min(zc)
  rho = np.empty((M+1, N+1))
  u = np.empty((M+1, N+1))
  v = np.empty((M+1, N+1))
  Pressure = np.empty((M+1, N+1)) 
  height = np.empty((M+1, N+1))
  engyDens = np.empty((M+1, N+1))
#  q = np.empty((3, N+2))
  s = np.empty((M+1, N+1)) 
  c_gas = np.empty((M+1, N+1))
  Ndens = np.empty((M+1, N+1))
  gamma1 = gamma - 1  



######_________Gas Dynamics_______________
  if (sim_type == 'gasDyno'):

    if (sim_details == 'shockTube'):
      print('====SHOCK TUBE=====')
      for m in range(M+1):
        for n in range(N+1):
          if (m <= (M-1)/2. and n <= (N-1)/2.): #top left quadrant
            u[m,n] = 1.206
            v[m,n] = 0.
            rho[m,n] = 0.5323
            Pressure[m,n] = 0.3 #N * k_B * Temp[i]
            engyDens[m,n] = (5./2.)* Pressure[m,n] + 0.5 * rho[m,n] *  u[m,n] * u[m,n]
            c_gas[m,n] = np.sqrt(gamma*Pressure[m,n]/rho[m,n])

          elif (m <= (M-1)/2. and n > (N-1)/2.): # top right quadrant
            u[m,n] = 1.206
            v[m,n] = 1.206
            rho[m,n] = 0.138
            Pressure[m,n] = 0.029 #N * k_B * Temp[i]
            engyDens[m,n] = (5./2.)* Pressure[m,n] + 0.5 * rho[m,n] *  u[m,n] * u[m,n]
            c_gas[m,n] = np.sqrt(gamma*Pressure[m,n]/rho[m,n])
   
          elif (m > (M-1)/2. and n <= (N-1)/2.): #bottom left quandrant
            u[m,n] = 0.
            v[m,n] = 0.
            rho[m,n] = 1.5
            Pressure[m,n] = 1.5 #N * k_B * Temp[i]
            engyDens[m,n] = (5./2.)* Pressure[m,n] + 0.5 * rho[m,n] *  u[m,n] * u[m,n]
            c_gas[m,n] = np.sqrt(gamma*Pressure[m,n]/rho[m,n])

          elif (m > (M-1)/2. and n > (N-1)/2.): #bottom right quandrant
            u[m,n] = 0.
            v[m,n] = 1.206
            rho[m,n] = 0.5323
            Pressure[m,n] = 0.3 #N * k_B * Temp[i]
            engyDens[m,n] = (5./2.)* Pressure[m,n] + 0.5 * rho[m,n] *  u[m,n] * u[m,n]
            c_gas[m,n] = np.sqrt(gamma*Pressure[m,n]/rho[m,n])


          q[0, m+1, n+1] = rho[m,n]
          q[1, m+1, n+1] = rho[m,n]*u[m,n] #u[j]
          q[2, m+1, n+1] = rho[m,n]*v[m,n]
          q[3, m+1, n+1] = engyDens[m,n] #Pressure[j]


    if (sim_details == 'constant'):
      print('====CONSTANT PRESSURE====')
          
      for j in range(N):
        if (j <= (N-1)/2.): #LEFT STATE
          u.append(0.)
          rho.append(1.)
          Pressure.append(1.) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] *  u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))
  
        elif (j > (N-1)/2.): #RIGHT STATE
          u.append(0.)
          rho.append(1.)
          Pressure.append(1.) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] * u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))

        q[0, j+1] = rho[j]
        q[1, j+1] = rho[j]*u[j] #u[j]
        q[2, j+1] = rho[j]*v[j]
        q[3, j+1] = engyDens[j] #Pressure[j]



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
        q[2, j+1] = rho[j]*v[j]
        q[3, j+1] = engyDens[j] #Pressure[j]

  return q, u, v, Pressure#, Ndens#,, Z, K, rho

def Mat_ghostCells(q, M, N, BC, DIM):

  if(DIM == 1):

    if(BC == 'extrap'):
      q[0,0] = q[0,1]
      q[0,-1] = q[0,-2]
      q[1,0] = q[1,1]
      q[1,-1] = q[1,-2]
      if q.shape[0] == 3:
        q[2,0] = q[2,1]
        q[2,-1] = q[2,-2]
  
    if(BC == 'wall'):
      q[0,0] = q[0,1]
      q[0,-1] = q[0,-2]
      q[1,0] = -q[1,1]
      q[1,-1] = q[1,-2]
      if q.shape[0] == 3:
        q[2,0] = q[2,1]
        q[2,-1] = q[2,-2]

  if(DIM == 2):

    if(BC == 'extrap'):
#---r spatial variable
      q[0,0,:] = q[0,1,:]
      q[0,-1,:] = q[0,-2,:]
      q[1,0,:] = q[1,1,:]
      q[1,-1,:] = q[1,-2,:]
      q[2,0,:] = q[2,1,:]
      q[2,-1,:] = q[2,-2,:]
      q[3,0,:] = q[3,1,:]
      q[3,-1,:] = q[3,-2,:]
#---z spatial variable
      q[0,:,0] = q[0,:,1]
      q[0,:,-1] = q[0,:,-2]
      q[1,:,0] = q[1,:,1]
      q[1,:,-1] = q[1,:,-2]
      q[2,:,0] = q[2,:,1]
      q[2,:,-1] = q[2,:,-2]
      q[3,:,0] = q[3,:,1]
      q[3,:,-1] = q[3,:,-2]
  
    if(BC == 'wall'):
  
      q[0,0] = q[0,1]
      q[0,-1] = q[0,-2]
      q[1,0] = -q[1,1]
      q[1,-1] = q[1,-2]
      if q.shape[0] == 3:
        q[2,0] = q[2,1]
        q[2,-1] = q[2,-2]

  return q

def fill_ghosts(q, N, Type):
#  for i in range(ng-1):
  if (Type == 'Wall'):
    #left wall
    q.insert(0, -q[1])
    #right wall
    q.insert(len(q), -q[-1])
  if (Type == 'Periodic'):
    #left Boundary
    q.insert(0, q[-1])
    #right Boundary 
    q.insert(len(q), q[1])
  
  if (Type == 'extrap'):
    q.insert(0, q[0])
    q.insert(len(q), q[-1])
  return q


def JumpSplit(q, sim_type, DIM, SPATv): 
  dq1 = 0.
  dq2 = 0.
  dq3 = 0.
  dq4 = 0.
  if (DIM == 1):
#  for i in range(1,N+1):
    dq1 = q[0,1:] - q[0,:-1] #q[0,i] - q[0,i-1]
    dq2 = q[1,1:] - q[1,:-1] #q[1,i] - q[1,i-1]
    if (sim_type == 'gasDyno'):   
      dq3 = q[2,1:] - q[2,:-1] #q[2,i] - q[2,i-1]
    else:
      dq3 = 0.

  if (DIM == 2):
    if (SPATv == 'r'):
      dq1 = q[0,1:,:] - q[0,:-1,:] 
      dq2 = q[1,1:,:] - q[1,:-1,:] 
      if (sim_type == 'gasDyno'):   
        dq3 = q[2,1:,:] - q[2,:-1,:]
        dq4 = q[3,1:,:] - q[3,:-1,:] 
    if (SPATv == 'z'):
      dq1 = q[0,:,1:] - q[0,:,:-1] 
      dq2 = q[1,:,1:] - q[1,:,:-1] 
      if (sim_type == 'gasDyno'):   
        dq3 = q[2,1:,:] - q[2,:-1,:]
        dq4 = q[3,1:,:] - q[3,:-1,:] 
  else:
    dq3 = 0.
    dq4 = 0.    

  return dq1,dq2,dq3,dq4
