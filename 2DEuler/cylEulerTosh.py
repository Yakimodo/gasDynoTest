###====This is a test euler program, if broken try eulerEqs.py

import numpy as np
import matplotlib.pyplot as plt
import twoD_function_list as fl

EFIX = 0
MOVIE = 0
HEAT = 0
RADIAL = 0
DIM = 2
CYL = 0
TITLE = 'shockTube'
#TITLE = 'csun-BalbasICs'
#TITLE = 'constant'
#TITLE = 'twoShock'


M = 100
N = 100
L = 1. #domain_length
CFL = 0.75 #Courant-Fredrichs-Lewy condition
dr = L / (2.*(M-1)) #spatial res in r
dz = L / (2.*(N-1)) #spatial res in z
dt = CFL*dr #time step
print('dt = ' + str(dt))
T_begin = 1 #step number [important in determining xi = x/t ----> (in this case) --> xc/(t*dt)]
T_end = 501 #number of steps
gamma = 1.4 #gravity
kB = 1.38e-23
rc = np.linspace(0., L, M+2)
zc = np.linspace(0., L, N+2)


eqNum = 4 #number of equations in system (Cyl Euler = 4) (Cyl Euler Isentropic = 3) (Euler Cartesian = 3) (SW = 2)
waveNum = eqNum #number of waves per Riemann Prob

q_zeros = np.empty((eqNum, M+2, N+2))
q_zeros[:,:,:] = 0.
W = np.zeros((eqNum, waveNum, M+1))
Wf = np.zeros((eqNum, waveNum, M+1, N+1))
Wg = np.zeros((eqNum, waveNum, M+1, N+1))
R = np.zeros((eqNum, waveNum, M+1))
Rf = np.zeros((eqNum, waveNum, M+1, N+1))
Rg = np.zeros((eqNum, waveNum, M+1, N+1))
u_hat = np.zeros((M+1, N+1))
v_hat = np.zeros((M+1, N+1))
H_hat = np.zeros((M+1, N+1))
H_hatR = np.zeros((M+1, N+1))
H_hatZ = np.zeros((M+1, N+1))
c_hat = np.zeros((M+1, N+1))
c_hatR = np.zeros((M+1, N+1))
c_hatZ = np.zeros((M+1, N+1))
a = np.zeros((eqNum, M+1))
af = np.zeros((eqNum, M+1, N+1))
ag = np.zeros((eqNum, M+1, N+1))
s = np.zeros((eqNum, M+1))
sf = np.zeros((eqNum, M+1, N+1))
sg = np.zeros((eqNum, M+1, N+1))
amdq = np.zeros((eqNum, M+1))
apdq = np.zeros((eqNum, M+1))
amdqf = np.zeros((eqNum, M+1, N+1))
apdqf = np.zeros((eqNum, M+1, N+1))
amdqg = np.zeros((eqNum, M+1, N+1))
apdqg = np.zeros((eqNum, M+1, N+1))
Q_heat = []

q, u, v, Pressure = fl.initialConditions('gasDyno', TITLE, rc, zc, M, N, q_zeros)

for t in range(T_begin, T_end+1):

  print('Step number = ' + str(t))
  q_new_a = np.zeros((q.shape[0],q.shape[1], q.shape[2]))
  q_new_b = np.zeros((q.shape[0],q.shape[1], q.shape[2]))
  q_new_c = np.zeros((q.shape[0],q.shape[1], q.shape[2]))

  q = fl.Mat_ghostCells(q, M, N, 'extrap', DIM)
#  print('q filled ghost cells = \n' + str(np.round(q,2)))

#  Q_heat = []
#  if (HEAT == 1):
#    for i in range(N+1):
#      Q_heat.append(25.*np.exp(-0.05*t)*EfOld[i]*JOld[i])


  rholR = q[0,1:-1,:-1] #(M, N+1) 
  rhorR = q[0,1:-1,1:]    #
  momlR = q[1,1:-1,:-1]   #
  momrR = q[1,1:-1,1:]    #
  ul = momlR/rholR        #
  ur = momrR/rhorR        #
  ElR = q[3,1:-1,:-1]     #
  ErR = q[3,1:-1,1:]  #(M, N+1)
  u = q[1,:,:]/q[0,:,:] #(M+2, N+2)

  rholZ = q[0,:-1,1:-1] #(M+1, N)
  rhorZ = q[0,1:,1:-1]   #
  momlZ = q[2,:-1,1:-1]  #
  momrZ = q[2,1:,1:-1]   #
  vl = momlZ/rholZ       #
  vr = momrZ/rhorZ       # 
  ElZ = q[3,:-1,1:-1]    #
  ErZ = q[3,1:,1:-1]    #(M+1, N)
  v = q[2,:,:]/q[0,:,:] #(M+2, N+2)

  P = (gamma-1)*(q[3,:,:]-0.5*q[0,:,:]*(u**2+v**2)) #(M+2, N+2)
  plR = P[1:-1,:-1] #(M,N+1)
  prR = P[1:-1,1:] 
  plZ = P[:-1,1:-1] #(M+1,N)
  prZ = P[1:,1:-1]
  rho = q[0,:,:] #(M+2,N+2)
  P[0,0] = 0.
  P[-1,0] = 0.
  P[0,-1] = 0.
  P[-1,-1] = 0.

#  if (MOVIE == 0):
#    R,Z = np.meshgrid(rc[1:-1], zc[1:-1])
#    plt.subplot(2,2,1)
#    plt.suptitle('time = ' + str(dt*t))
#    plt.title('mass density' )
#    cp = plt.contourf(rho[:,:])
#    plt.colorbar(cp)
#    plt.subplot(2,2,2)
#    plt.title('Pressure' )
#    cp = plt.contourf(P[:,:])
#    plt.colorbar(cp)
#    plt.show()

###-----2D Cylindrical-----------------------
  for m in range(M+1): #r
    for n in range(N): #z
      v_hat[m,n] = (np.sqrt(rholZ[m,n])*vl[m,n] + np.sqrt(rhorZ[m,n])*vr[m,n]) / (np.sqrt(rholZ[m,n]) + np.sqrt(rhorZ[m,n])) #(M+1,N)

#      H_hatZ[m,n] = ((ElZ[m,n] + plZ[m,n]) / np.sqrt(rholZ[m,n]) + (ErZ[m,n] + prZ[m,n]) / np.sqrt(rhorZ[m,n])) / (np.sqrt(rholZ[m,n]) + np.sqrt(rhorZ[m,n])) 
#      c_hatZ[m,n] = np.sqrt((gamma-1)*(H_hatZ[m,n]-0.5*(u_hat[m,n]**2+v_hat[m,n]**2)))

  for m in range(M): #r
    for n in range(N+1): #z
      u_hat[m,n] = (np.sqrt(rholR[m,n])*ul[m,n] + np.sqrt(rhorR[m,n])*ur[m,n]) / (np.sqrt(rholR[m,n]) + np.sqrt(rhorR[m,n])) #(M+1,N+1)
      H_hatR[m,n] = ((ElR[m,n] + plR[m,n]) / np.sqrt(rholR[m,n]) + (ErR[m,n] + prR[m,n]) / np.sqrt(rhorR[m,n])) / (np.sqrt(rholR[m,n]) + np.sqrt(rhorR[m,n])) #(M+1,N+1)
      c_hatR[m,n] = np.sqrt((gamma-1)*(H_hatR[m,n]-0.5*(u_hat[m,n]**2 + v_hat[m,n]**2))) #(M+1,N+1)

  
  dq1R, dq2R, dq3R, dq4R = fl.JumpSplit(q, 'gasDyno', DIM, 'r')

  for m in range(M): #r
    for n in range(N+1): #z
###-----r subroutine
      af[2,m,n] = dq3R[m,n] - dq1R[m,n] * v_hat[m,n]
      af[1,m,n] = ((gamma-1)/c_hatR[m,n]**2) * (dq1R[m,n]*(H_hatR[m,n]-u_hat[m,n]**2) - u_hat[m,n]*dq2R[m,n] - dq4R[m,n] -  (dq3R[m,n] - v_hat[m,n]*dq1R[m,n])*v_hat[m,n])
      af[0,m,n] = 1/(2.*c_hatR[m,n]) * (dq1R[m,n]*(u_hat[m,n]+c_hatR[m,n]) - dq2R[m,n] - c_hatR[m,n]*af[1,m,n])
      af[3,m,n] = dq1R[m,n] - af[0,m,n] - af[1,m,n]
   
      sf[0,m,n] = u_hat[m,n] - c_hatR[m,n]
      sf[1,m,n] = u_hat[m,n]
      sf[2,m,n] = u_hat[m,n]
      sf[3,m,n] = u_hat[m,n] + c_hatR[m,n]
    
      Rf[0,0,m,n] = 1.
      Rf[0,1,m,n] = u_hat[m,n] - c_hatR[m,n]
      Rf[0,2,m,n] = v_hat[m,n]
      Rf[0,3,m,n] = H_hatR[m,n] - u_hat[m,n]*c_hatR[m,n]

      Rf[1,0,m,n] = 1.
      Rf[1,1,m,n] = u_hat[m,n]
      Rf[1,2,m,n] = v_hat[m,n]
      Rf[1,3,m,n] = 0.5 *( u_hat[m,n]**2 + v_hat[m,n]**2 )

      Rf[2,0,m,n] = 0.
      Rf[2,1,m,n] = 0.
      Rf[2,2,m,n] = 1.
      Rf[2,3,m,n] = v_hat[m,n]

      Rf[3,0,m,n] = 1.
      Rf[3,1,m,n] = u_hat[m,n] + c_hatR[m,n]
      Rf[3,2,m,n] = v_hat[m,n]
      Rf[3,3,m,n] = H_hatR[m,n] + u_hat[m,n]*c_hatR[m,n]


  Wf[:,:,:,:] = 0.
  for j in range(eqNum):
    for w in range(waveNum):
      for m in range(M+1): #spatial variable r
        for n in range(N+1): #spatial variable z
          Wf[w,j,m,n] = af[w,m,n] * Rf[w,j,m,n]


  amdqf[:,:,:] = 0.
  apdqf[:,:,:] = 0.

  for j in range(eqNum): 
    for m in range(M):
      for n in range(N):
        for w in range(waveNum):
          amdqf[j,m,n] += min(sf[w,m,n+1],0) * Wf[w,j,m,n+1]
          apdqf[j,m,n] += max(sf[w,m,n],0) * Wf[w,j,m,n]   

#  print('amdqf = \n' + str(np.round(amdqf,2))) 
#  print('apdqf = \n' + str(np.round(apdqf,2))) 

  q_new_a[:,:,:] = 0.
  for j in range(eqNum): 
    for m in range(M): #q = q.shape[0],q.shape[1],q.shape[2] (eqnum,n+2,m+2)
      for n in range(N+1):
        q_new_a[j,m+1,n+1] = q[j,m+1,n+1] - (1./2.) * dt/dr * (amdqf[j,m,n] + apdqf[j,m,n])   
        

#  print('prior to filling ghost  q_new_a = \n'  + str(q_new_a))
  q_new_a = fl.Mat_ghostCells(q_new_a, M, N, 'extrap', DIM)
#  print('q_new_a = \n' + str(np.round(q_new_a,2)))

  rholR = q_new_a[0,1:-1,:-1] #(M, N+1) 
  rhorR = q_new_a[0,1:-1,1:]    #
  momlR = q_new_a[1,1:-1,:-1]   #
  momrR = q_new_a[1,1:-1,1:]    #
  ul = momlR/rholR        #
  ur = momrR/rhorR        #
  ElR = q_new_a[3,1:-1,:-1]     #
  ErR = q_new_a[3,1:-1,1:]  #(M, N+1)
  u = q_new_a[1,:,:]/q_new_a[0,:,:] #(M+2, N+2)

  rholZ = q_new_a[0,:-1,1:-1] #(M+1, N)
  rhorZ = q_new_a[0,1:,1:-1]   #
  momlZ = q_new_a[2,:-1,1:-1]  #
  momrZ = q_new_a[2,1:,1:-1]   #
  vl = momlZ/rholZ       #
  vr = momrZ/rhorZ       # 
  ElZ = q_new_a[3,:-1,1:-1]    #
  ErZ = q_new_a[3,1:,1:-1]    #(M+1, N)
  v = q_new_a[2,:,:]/q_new_a[0,:,:] #(M+2, N+2)

  P = (gamma-1)*(q_new_a[3,:,:]-0.5*q_new_a[0,:,:]*(u**2+v**2)) #(M+2, N+2)
  plR = P[1:-1,:-1] #(M,N+1)
  prR = P[1:-1,1:] 
  plZ = P[:-1,1:-1] #(M+1,N)
  prZ = P[1:,1:-1]
  rho = q_new_a[0,:,:] #(M+2,N+2)
  P[0,0] = 0.
  P[-1,0] = 0.
  P[0,-1] = 0.
  P[-1,-1] = 0.


#  if (MOVIE == 0):
#    plt.subplot(2,2,1)
#    plt.suptitle('r subroutine time = ' + str(dt*t))
#    plt.title('mass density' )
#    cp = plt.contourf(rho[:,:])
#    plt.colorbar(cp)
#    plt.subplot(2,2,2)
#    plt.title('Pressure' )
#    cp = plt.contourf(P[:,:])
#    plt.colorbar(cp)
#    plt.show()

  dq1Z, dq2Z, dq3Z, dq4Z = fl.JumpSplit(q_new_a, 'gasDyno', DIM, 'z')

  for m in range(M): #r
    for n in range(N+1): #z
      u_hat[m,n] = (np.sqrt(rholR[m,n])*ul[m,n] + np.sqrt(rhorR[m,n])*ur[m,n]) / (np.sqrt(rholR[m,n]) + np.sqrt(rhorR[m,n]))
#      H_hatR[m,n] = ((ElR[m,n] + plR[m,n]) / np.sqrt(rholR[m,n]) + (ErR[m,n] + prR[m,n]) / np.sqrt(rhorR[m,n])) / (np.sqrt(rholR[m,n]) + np.sqrt(rhorR[m,n])) 
#      c_hatR[m,n] = np.sqrt((gamma-1)*(H_hatR[m,n]-0.5*(u_hat[m,n]**2+v_hat[m,n]**2)))

  for m in range(M+1): #r
    for n in range(N): #z
#      if (DIM == 2): 
      v_hat[m,n] = (np.sqrt(rholZ[m,n])*vl[m,n] + np.sqrt(rhorZ[m,n])*vr[m,n]) / (np.sqrt(rholZ[m,n]) + np.sqrt(rhorZ[m,n]))
      H_hatZ[m,n] = ((ElZ[m,n] + plZ[m,n]) / np.sqrt(rholZ[m,n]) + (ErZ[m,n] + prZ[m,n]) / np.sqrt(rhorZ[m,n])) / (np.sqrt(rholZ[m,n]) + np.sqrt(rhorZ[m,n])) 
      c_hatZ[m,n] = np.sqrt((gamma-1)*(H_hatZ[m,n]-0.5*(u_hat[m,n]**2+v_hat[m,n]**2)))

###----z subroutine
  for m in range(M+1): #r
    for n in range(N): #z
      ag[2,m,n] = dq2Z[m,n] - dq1Z[m,n]*u_hat[m,n]
      ag[1,m,n] = ((gamma-1)/c_hatZ[m,n]**2) * (dq1Z[m,n]*(H_hatZ[m,n]-v_hat[m,n]**2) - v_hat[m,n]*dq3Z[m,n] - dq4Z[m,n] - (dq2Z[m,n] - u_hat[m,n]*dq1Z[m,n])*u_hat[m,n])
      ag[0,m,n] = 1/(2.*c_hatZ[m,n]) * (dq1Z[m,n]* (v_hat[m,n]+c_hatZ[m,n]) - dq3Z[m,n] - c_hatZ[m,n]*ag[1,m,n])
      ag[3,m,n] = dq1Z[m,n] - ag[0,m,n] - ag[1,m,n]

      sg[0,m,n] = v_hat[m,n] - c_hatZ[m,n]
      sg[1,m,n] = v_hat[m,n]
      sg[2,m,n] = v_hat[m,n]
      sg[3,m,n] = v_hat[m,n] + c_hatZ[m,n]
    
      Rg[0,0,m,n] = 1.
      Rg[0,1,m,n] = u_hat[m,n]
      Rg[0,2,m,n] = v_hat[m,n] - c_hatZ[m,n]
      Rg[0,3,m,n] = H_hatZ[m,n] - v_hat[m,n]*c_hatZ[m,n]

      Rg[1,0,m,n] = 1.
      Rg[1,1,m,n] = u_hat[m,n]
      Rg[1,2,m,n] = v_hat[m,n]
      Rg[1,3,m,n] = 0.5 *( u_hat[m,n]**2 + v_hat[m,n]**2 )

      Rg[2,0,m,n] = 0.
      Rg[2,1,m,n] = 1.
      Rg[2,2,m,n] = 0.
      Rg[2,3,m,n] = u_hat[m,n]

      Rg[3,0,m,n] = 1.
      Rg[3,1,m,n] = u_hat[m,n]
      Rg[3,2,m,n] = v_hat[m,n] + c_hatZ[m,n]
      Rg[3,3,m,n] = H_hatZ[m,n] + u_hat[m,n]*c_hatZ[m,n]

#  print('ag = ' + str(ag))
#  print('Rg = ' + str(Rg))

  Wg[:,:,:,:] = 0.
  for j in range(eqNum):
    for w in range(waveNum):
      for m in range(M+1): #spatial variable r
        for n in range(N+1): #spatial variable z
          Wg[w,j,m,n] = ag[w,m,n] * Rg[w,j,m,n]

  amdqg[:,:,:] = 0.
  apdqg[:,:,:] = 0.
 
#      print('sg = ' + str(sg))
#      print('Wg = ' + str(Wg))
   
  for j in range(eqNum): 
    for m in range(M):
      for n in range(N):
        for w in range(waveNum):
          amdqg[j,m,n] += min(sg[w,m+1,n],0) * Wg[w,j,m+1,n]
          apdqg[j,m,n] += max(sg[w,m,n],0) * Wg[w,j,m,n]   
 
#  print('amdqg = \n' + str(np.round(amdqg,2))) 
#  print('apdqg = \n' + str(np.round(apdqg,2))) 

  q_new_b[:,:,:] = 0.
  for j in range(eqNum): 
    for m in range(M+1): #q = q.shape[0],q.shape[1],q.shape[2] (eqnum,n+2,m+2)
      for n in range(N):
        q_new_b[j,m+1,n+1] = q_new_a[j,m+1,n+1] - (1./2.) * dt/dz * (amdqg[j,m,n] + apdqg[j,m,n])  
    
#  print('q_new_b prior to ghosts = \n' + str(np.round(q_new_b,2)))
                ##Add ghost cell filling function
  q_new_b = fl.Mat_ghostCells(q_new_b, M, N, 'extrap', DIM)
#  print('q_new_b = \n' + str(np.round(q_new_b,2)))
#  print('q_new_a - q_new_b = ' + str(q_new_a - q_new_b))
  if (CYL == 1):
    q_new_c = q_new_b
    for j in range(eqNum): 
      for m in range(m): #q = q.shape[0],q.shape[1],q.shape[2] (eqnum,n+2,m+2)
        for n in range(n):
          if (j == 0):
            q_new_c[j,m+1,n+1] = q_new_b[j,m,n+1] - (1./3.)*dt * (1/rc[m+1] * 2.*q_new_b[1,m,n+1])
          elif (j == 1):
            q_new_c[j,m+1,n+1] = q_new_b[j,m, n+1] - (1./3.)*dt * (1/rc[m+1] * (2.*q_new_b[1,m,n+1]**2/q_new_b[0,m,n+1]))
          elif (j == 2):
            q_new_c[j,m+1,n+1] = q_new_b[j,m, n+1] - (1./3.)*dt * (1/rc[m+1] * (2.*q_new_b[2,m,n+1]*q_new_b[1,m,n+1]/q_new_b[0,m,n+1]))
          elif (j == 3):
            q_new_c[j,m+1,n+1] = q_new_b[j,m, n+1] - (1./3.)*dt * (1/rc[m+1] * ((q_new_b[3,m,n+1]+P[m+1,n+1])*(q_new_b[1,m,n+1]/q_new_b[0,n+1])))
        
 
  q_new_b = fl.Mat_ghostCells(q_new_b, M, N, 'extrap', DIM)
#  print('q_new_b = \n' + str(q_new_b))
#  print('q - q_new_b = ' + str(q-q_new_b))

  rholR = q_new_b[0,1:-1,:-1] #(M, N+1) 
  rhorR = q_new_b[0,1:-1,1:]    #
  momlR = q_new_b[1,1:-1,:-1]   #
  momrR = q_new_b[1,1:-1,1:]    #
  ul = momlR/rholR        #
  ur = momrR/rhorR        #
  ElR = q_new_b[3,1:-1,:-1]     #
  ErR = q_new_b[3,1:-1,1:]  #(M, N+1)
  u = q_new_b[1,:,:]/q_new_b[0,:,:] #(M+2, N+2)

  rholZ = q_new_b[0,:-1,1:-1] #(M+1, N)
  rhorZ = q_new_b[0,1:,1:-1]   #
  momlZ = q_new_b[2,:-1,1:-1]  #
  momrZ = q_new_b[2,1:,1:-1]   #
  vl = momlZ/rholZ       #
  vr = momrZ/rhorZ       # 
  ElZ = q_new_b[3,:-1,1:-1]    #
  ErZ = q_new_b[3,1:,1:-1]    #(M+1, N)
  v = q_new_b[2,:,:]/q_new_b[0,:,:] #(M+2, N+2)

  P = (gamma-1)*(q_new_b[3,:,:]-0.5*q_new_b[0,:,:]*(u**2+v**2)) #(M+2, N+2)
  plR = P[1:-1,:-1] #(M,N+1)
  prR = P[1:-1,1:] 
  plZ = P[:-1,1:-1] #(M+1,N)
  prZ = P[1:,1:-1]
  rho = q_new_b[0,:,:] #(M+2,N+2)
  P[0,0] = 0.
  P[-1,0] = 0.
  P[0,-1] = 0.
  P[-1,-1] = 0.

###___PLOTTING___###

  if(MOVIE == 0):  
#    if(t%75 == 0 or t == 1): 
    if(t%25 == 0 or t == 1): 
      plt.subplot(2,2,1)
      plt.suptitle('time = ' + str(dt*t))
      plt.title('mass density' )
      cp = plt.contourf(rho[:,:])
      plt.colorbar(cp)
      plt.subplot(2,2,2)
      plt.title('Pressure' )
      cp = plt.contourf(P[:,:])
      plt.colorbar(cp)
      plt.show()

#      plt.figure(1)
#      R,Z = np.meshgrid(rc, zc)
#      plt.title(' time  = ' + str(dt*t) )
#      cp = plt.contourf(rho, label = 'density')
#      plt.colorbar(cp)
#      plt.axvline(x=0., color='k', linestyle='--')  #(max(rc)/2.))
#      plt.legend()
#      plt.show()

  if(MOVIE == 1):  
    plt.clf()
    plt.subplot(2,2,1)
    plt.suptitle('z subroutine time = ' + str(dt*t))
    plt.title('mass density' )
    cp = plt.contourf(rho[:,:])
    plt.colorbar(cp)
    plt.subplot(2,2,2)
    plt.title('Pressure' )
    cp = plt.contourf(P[:,:])
    plt.colorbar(cp)



#    if(TITLE == 'constant'):
#      plt.axis([rc.min(),rc.max(),-1.1,2.1])
    plt.pause(1.5)

  if (HEAT == 0):
    if (CYL == 0):
      q = q_new_b
    elif (CYL == 1):
      q = q_new_c








