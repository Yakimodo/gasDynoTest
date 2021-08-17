import numpy as np
import matplotlib.pyplot as plt

kB = 1.38e-23
gamma = 1.4 #air == 7/5
beta = 0.2

def initialConditions(sim_type, sim_details, xc, N, q, DIM):
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
  Ndens = []
  Z = []
  gamma1 = gamma - 1  
  g = 1. # gravity for shallow water eqs.



######_________Gas Dynamics_______________
  if (sim_type == 'gasDyno'):

    if (sim_details == 'shockTube'):
      print('====SHOCK TUBE=====')
      for j in range(N):
        if (j <= (N-1)/2.): #LEFT STATE
          u.append(0.)
          v.append(0.)
          rho.append(3.)
          Pressure.append(3.) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] *  u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))
#          Ndens.append(0.01*np.exp(-np.abs((xc[j]-xc[np.int((N-1)/2)]))**2 / 0.2))
 
        elif (j > (N-1)/2.): #RIGHT STATE
          u.append(0.)
          v.append(0.)
          rho.append(1.)
          Pressure.append(1.) #N * k_B * Temp[i])
          engyDens.append((5./2.)* Pressure[j] + 0.5 * rho[j] * u[j] * u[j])
          c_gas.append(np.sqrt(gamma*Pressure[j]/rho[j]))
#          Ndens.append(0.01**np.exp(-np.abs((xc[j]-xc[np.int((N-1)/2)]))**2 / 0.2))

        q[0, j+1] = rho[j]
        q[1, j+1] = rho[j]*u[j] #u[j]
        q[2, j+1] = rho[j]*v[j]
        q[3, j+1] = engyDens[j] #Pressure[j]


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

def Mat_ghostCells(q, N, BC):

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


def JumpSplit(q, N, sim_type, DIM): 
  dq1 = 0.
  dq2 = 0.
  dq3 = 0.
  dq4 = 0.
#  for i in range(1,N+1):
  dq1 = q[0,1:] - q[0,:-1] #q[0,i] - q[0,i-1]
  dq2 = q[1,1:] - q[1,:-1] #q[1,i] - q[1,i-1]
  if (sim_type == 'gasDyno'):   
    dq3 = q[2,1:] - q[2,:-1] #q[2,i] - q[2,i-1]
  if (DIM = 2):
    dq4 = q[3,1:] - q[3,:-1]   
  else:
    dq3 = 0.
    dq4 = 0.

  return dq1,dq2,dq3,dq4
