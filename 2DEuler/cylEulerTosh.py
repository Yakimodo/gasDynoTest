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
#TITLE = 'constant'
#TITLE = 'twoShock'


M = 10
N = 10
L = 1. #domain_length
CFL = 0.6 #Courant-Fredrichs-Lewy condition
dr = L / (2.*(M-1)) #spatial res in r
dz = L / (2.*(N-1)) #spatial res in z
dt = CFL*dr #time step
T_begin = 1 #step number [important in determining xi = x/t ----> (in this case) --> xc/(t*dt)]
T_end = 500 #number of steps
gamma = 1.4 #gravity
kB = 1.38e-23
rc = np.linspace(0., L, M+1)
zc = np.linspace(0., L, N+1)


eqNum = 4 #number of equations in system (Cyl Euler = 4) (Cyl Euler Isentropic = 3) (Euler Cartesian = 3) (SW = 2)
waveNum = eqNum #number of waves per Riemann Prob

q_empty = np.empty((eqNum, M+2, N+2))
W = np.empty((eqNum, waveNum, M+1))
Wf = np.empty((eqNum, waveNum, M+1, N+1))
Wg = np.empty((eqNum, waveNum, M+1, N+1))
R = np.empty((eqNum, waveNum, M+1))
Rf = np.empty((eqNum, waveNum, M+1, N+1))
Rg = np.empty((eqNum, waveNum, M+1, N+1))
u_hat = np.empty((M+1, N+1))
v_hat = np.empty((M+1, N+1))
H_hat = np.empty((M+1, N+1))
H_hatR = np.empty((M+1, N+1))
H_hatZ = np.empty((M+1, N+1))
c_hat = np.empty((M+1, N+1))
c_hatR = np.empty((M+1, N+1))
c_hatZ = np.empty((M+1, N+1))
a = np.empty((eqNum, M+1))
af = np.empty((eqNum, M+1, N+1))
ag = np.empty((eqNum, M+1, N+1))
s = np.empty((eqNum, M+1))
sf = np.empty((eqNum, M+1, N+1))
sg = np.empty((eqNum, M+1, N+1))
amdq = np.empty((eqNum, M+1))
apdq = np.empty((eqNum, M+1))
amdqf = np.empty((eqNum, M+1, N+1))
apdqf = np.empty((eqNum, M+1, N+1))
amdqg = np.empty((eqNum, M+1, N+1))
apdqg = np.empty((eqNum, M+1, N+1))
Q_heat = []

q, u, v, P = fl.initialConditions('gasDyno', TITLE, rc, zc, M, N, q_empty)
#q, u, P = fl.initialConditions('gasDyno', 'constant', rc, N, q_empty)
#q, u, P = fl.initialConditions('gasDyno', 'twoShock', rc, N, q_empty)
#  q0 = density
#  q1 = velocity
#  q2 = pressure

#Ndens = []
#for j in range(N):
#  Ndens.append(1.e20*np.exp(-np.abs((rc[j]-rc[np.int((N-1)/2)]))**2 / 0.2))
#
#Ndens = fl.fill_ghosts(Ndens, N, 'extrap')
#Ndens = np.array(Ndens)


for t in range(T_begin, T_end+1):

  q_new_a = np.empty((q.shape[0],q.shape[1], q.shape[2]))
  q_new_b = np.empty((q.shape[0],q.shape[1], q.shape[2]))
  q_new_c = np.empty((q.shape[0],q.shape[1], q.shape[2]))

#  if (RADIAL == 0):
  q = fl.Mat_ghostCells(q, M, N, 'extrap', DIM)

  Q_heat = []
  if (HEAT == 1):
    for i in range(N+1):
      Q_heat.append(25.*np.exp(-0.05*t)*EfOld[i]*JOld[i])


  if (DIM == 2):
    rholR = q[0,:,:-1]
    rhorR = q[0,:,1:]
    rholZ = q[0,:-1,:]
    rhorZ = q[0,1:,:]
    momlR = q[1,:,:-1] #rhol * ul
    momrR = q[1,:,1:] #rhor * ur
    momlZ = q[2,:-1,:] #rhol * vl
    momrZ = q[2,1:,:] #rhor * vr
    ul = momlR/rholR #u[1,:-1]
    ur = momrR/rhorR #u[1,1:]
    vl = momlZ/rholZ #u[1,:-1]
    vr = momrZ/rhorZ #u[1,1:]
    u = q[1,:,:]/q[0,:,:]
    v = q[2,:,:]/q[0,:,:]
#    print(u.shape,v.shape)
    ElR = q[3,:,:-1]#pl/(gamma - 1) + 0.5 * rhol * ul**2 #q[1,:-1] / q[0,:-1]
    ErR = q[3,:,1:]#pr/(gamma - 1) + 0.5 * rhor * ur**2 #q[1,1:] / q[0,1:]  
    ElZ = q[3,:-1,:]#pl/(gamma - 1) + 0.5 * rhol * ul**2 #q[1,:-1] / q[0,:-1]
    ErZ = q[3,1:,:]#pr/(gamma - 1) + 0.5 * rhor * ur**2 #q[1,1:] / q[0,1:]  
#    print(ElR.shape,rholR.shape,ul.shape,vl.shape)
#    print(ErR.shape,rhorR.shape,ur.shape,vr.shape)
#    print('ElR = ' + str(ElR))
#    plR = (gamma-1)*(ElR[:,:-1]-0.5*rholR[:,:-1]*(ul[:,:-1]**2+vl[:-1,:]**2))#P[2,:-1]
#    prR = (gamma-1)*(ErR[:,:-1]-0.5*rhorR[:,:-1]*(ur[:,:-1]**2+vr[:-1,:]**2))#P[2,:-1]
#    plZ = (gamma-1)*(ElZ[:-1,:]-0.5*rholZ[:-1,:]*(ul[:,:-1]**2+vl[:-1,:]**2))#P[2,:-1]
#    prZ = (gamma-1)*(ErZ[:-1,:]-0.5*rhorZ[:-1,:]*(ur[:,:-1]**2+vr[:-1,:]**2))#P[2,:-1]
  #  Temp = 2./(5.*Ndens[:]*kB) * (q[2,:] - 0.5 * q[1,:]**2/q[0,:])
    P = (gamma-1)*(q[3,:,:]-0.5*q[0,:,:]*(u**2+v**2))
    plR = P[:,:-1]
    prR = P[:,1:]
    plZ = P[:-1,:]
    prZ = P[1:,:]
    rho = q[0,:,:]

#    print('plR = ' + str(plR))
#    print('rholR = ' + str(rholR))

#    clR = np.sqrt(gamma*plR/rholR)
#    crR = np.sqrt(gamma*prR/rhorR)
#    clZ = np.sqrt(gamma*plZ/rholZ)
#    crZ = np.sqrt(gamma*prZ/rhorZ)


###-----2D Cylindrical-----------------------
  for m in range(M+1): #r
    for n in range(N+1): #z

#      if (DIM == 2): 
      u_hat[m,n] = (np.sqrt(rholR[m,n])*ul[m,n] + np.sqrt(rhorR[m,n])*ur[m,n]) / (np.sqrt(rholR[m,n]) + np.sqrt(rhorR[m,n]))
      v_hat[m,n] = (np.sqrt(rholZ[m,n])*vl[m,n] + np.sqrt(rhorZ[m,n])*vr[m,n]) / (np.sqrt(rholZ[m,n]) + np.sqrt(rhorZ[m,n]))
      H_hatR[m,n] = ((ElR[m,n] + plR[m,n]) / np.sqrt(rholR[m,n]) + (ErR[m,n] + prR[m,n]) / np.sqrt(rhorR[m,n])) / (np.sqrt(rholR[m,n]) + np.sqrt(rhorR[m,n])) 
      H_hatZ[m,n] = ((ElZ[m,n] + plZ[m,n]) / np.sqrt(rholZ[m,n]) + (ErZ[m,n] + prZ[m,n]) / np.sqrt(rhorZ[m,n])) / (np.sqrt(rholZ[m,n]) + np.sqrt(rhorZ[m,n])) 
      c_hatR[m,n] = np.sqrt((gamma-1)*(H_hatR[m,n]-0.5*(u_hat[m,n]**2+v_hat[m,n]**2)))
      c_hatZ[m,n] = np.sqrt((gamma-1)*(H_hatZ[m,n]-0.5*(u_hat[m,n]**2+v_hat[m,n]**2)))

#  print(u_hat.shape,v_hat.shape,H_hatR.shape,H_hatZ.shape,c_hatR.shape,c_hatZ.shape)

  dq1R, dq2R, dq3R, dq4R = fl.JumpSplit(q, 'gasDyno', DIM, 'r')
#  print('dq1R = ' + str(dq1R) + '\ndq2R = ' + str(dq2R) + '\ndq3R = ' + str(dq3R) + '\ndq4R = ' + str(dq4R))
#  print('dq1Z = ' + str(dq1Z) + '\ndq2Z = ' + str(dq2Z) + '\ndq3Z = ' + str(dq3Z) + '\ndq4Z = ' + str(dq4Z))


  for m in range(M+1): #r
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

#  print('af = ' + str(af))
#  print('Rf = ' + str(Rf))
#  print('af[0,0,5] = ' + str(af[0,0,5]))
#  print('Rf[0,0,0,5] = ' + str(Rf[0,0,0,5]))

  Wf[:,:,:,:] = 0.
  for j in range(eqNum):
    for w in range(waveNum):
      for m in range(M+1): #spatial variable r
        for n in range(N+1): #spatial variable z
          Wf[w,j,m,n] = af[w,m,n] * Rf[w,j,m,n]
#          if (w == 0 and j == 0 and m == 0 and n == 5):
#            print('inside for loop  = ' + str(Wf[w,j,m,n]) + str(af[w,m,n]) + str(Rf[w,j,m,n]))

#  print('sf = ' + str(sf))
#  print('Wf[0,0,0,5] = ' + str(Wf[0,0,0,5]) + str(af[0,0,5]*Rf[0,0,0,5]))

  amdqf[:,:,:] = 0.
  apdqf[:,:,:] = 0.

  for j in range(eqNum): 
    for m in range(M-1):
      for n in range(N-1):
        for w in range(waveNum):
          amdqf[j,m,n] += min(sf[w,m+1,n],0) * Wf[w,j,m+1,n]
          apdqf[j,m,n] += max(sf[w,m,n],0) * Wf[w,j,m,n]   

  print('amdqf = ' + str(amdqf))
  print('apdqf = ' + str(apdqf))

  for j in range(eqNum): 
    for m in range(M+1): #q = q.shape[0],q.shape[1],q.shape[2] (eqnum,n+2,m+2)
      for n in range(N+1):
        q_new_a[j,m+1,n] = q[j,m+1,n] - (1./2.) * dt/dr * (amdqf[j,m,n] + apdqf[j,m,n])   


  plt.figure(2)
  R,Z = np.meshgrid(rc, zc)
  plt.title(' time  = ' + str(dt*t) )
#      plt.contourf(R, Z, 
  print('Pressure = ' + str(P))
#      print('rho = ' + str(rho))
#      print('u = ' + str(u))
#      print('v = ' + str(v))
  cp = plt.contourf(rho, label = 'density')
  plt.colorbar(cp)
  plt.show()
#  print('q = ' + str(q))
#  print('prior to filling ghost sq_new_a = ' + str(q_new_a))
  q_new_a = fl.Mat_ghostCells(q_new_a, M, N, 'extrap', DIM)
  print('q = ' + str(q))
  print('q_new_a = ' + str(q_new_a))
  print('q - q_new_a = ' + str(q-q_new_a))

  if (DIM == 2):
    rholR = q_new_a[0,:,:-1]
    rhorR = q_new_a[0,:,1:]
    rholZ = q_new_a[0,:-1,:]
    rhorZ = q_new_a[0,1:,:]
    print(rholR.shape,rhorR.shape,rholZ.shape,rhorZ.shape)
    momlR = q_new_a[1,:,:-1] #rhol * ul
    momrR = q_new_a[1,:,1:] #rhor * ur
    momlZ = q_new_a[2,:-1,:] #rhol * vl
    momrZ = q_new_a[2,1:,:] #rhor * vr
    ul = momlR/rholR #u[1,:-1]
    ur = momrR/rhorR #u[1,1:]
    vl = momlZ/rholZ #u[1,:-1]
    vr = momrZ/rhorZ #u[1,1:]
    u = q_new_a[1,:,:]/q_new_a[0,:,:]
    v = q_new_a[2,:,:]/q_new_a[0,:,:]
#    print(u.shape,v.shape)
    ElR = q_new_a[3,:,:-1]#pl/(gamma - 1) + 0.5 * rhol * ul**2 #q_new_a[1,:-1] / q_new_a[0,:-1]
    ErR = q_new_a[3,:,1:]#pr/(gamma - 1) + 0.5 * rhor * ur**2 #q_new_a[1,1:] / q_new_a[0,1:]  
    ElZ = q_new_a[3,:-1,:]#pl/(gamma - 1) + 0.5 * rhol * ul**2 #q_new_a[1,:-1] / q_new_a[0,:-1]
    ErZ = q_new_a[3,1:,:]#pr/(gamma - 1) + 0.5 * rhor * ur**2 #q_new_a[1,1:] / q_new_a[0,1:]  #    print(ElR.shape,rholR.shape,ul.shape,vl.shape)
#    print(ErR.shape,rhorR.shape,ur.shape,vr.shape)
#    print('ElR = ' + str(ElR))
#    plR = (gamma-1)*(ElR[:,:-1]-0.5*rholR[:,:-1]*(ul[:,:-1]**2+vl[:-1,:]**2))#P[2,:-1]
#    prR = (gamma-1)*(ErR[:,:-1]-0.5*rhorR[:,:-1]*(ur[:,:-1]**2+vr[:-1,:]**2))#P[2,:-1]
#    plZ = (gamma-1)*(ElZ[:-1,:]-0.5*rholZ[:-1,:]*(ul[:,:-1]**2+vl[:-1,:]**2))#P[2,:-1]
#    prZ = (gamma-1)*(ErZ[:-1,:]-0.5*rhorZ[:-1,:]*(ur[:,:-1]**2+vr[:-1,:]**2))#P[2,:-1]
  #  Temp = 2./(5.*Ndens[:]*kB) * (q[2,:] - 0.5 * q[1,:]**2/q[0,:])
    P = (gamma-1)*(q_new_a[3,:,:]-0.5*q_new_a[0,:,:]*(u**2+v**2))
    plR = P[:,:-1]
    prR = P[:,1:]
    plZ = P[:-1,:]
    prZ = P[1:,:]
    rho = q_new_a[0,:,:]

  dq1Z, dq2Z, dq3Z, dq4Z = fl.JumpSplit(q_new_a, 'gasDyno', DIM, 'z')

  for m in range(M+1): #r
    for n in range(N+1): #z

#      if (DIM == 2): 
      u_hat[m,n] = (np.sqrt(rholR[m,n])*ul[m,n] + np.sqrt(rhorR[m,n])*ur[m,n]) / (np.sqrt(rholR[m,n]) + np.sqrt(rhorR[m,n]))
      v_hat[m,n] = (np.sqrt(rholZ[m,n])*vl[m,n] + np.sqrt(rhorZ[m,n])*vr[m,n]) / (np.sqrt(rholZ[m,n]) + np.sqrt(rhorZ[m,n]))
      H_hatR[m,n] = ((ElR[m,n] + plR[m,n]) / np.sqrt(rholR[m,n]) + (ErR[m,n] + prR[m,n]) / np.sqrt(rhorR[m,n])) / (np.sqrt(rholR[m,n]) + np.sqrt(rhorR[m,n])) 
      H_hatZ[m,n] = ((ElZ[m,n] + plZ[m,n]) / np.sqrt(rholZ[m,n]) + (ErZ[m,n] + prZ[m,n]) / np.sqrt(rhorZ[m,n])) / (np.sqrt(rholZ[m,n]) + np.sqrt(rhorZ[m,n])) 
      c_hatR[m,n] = np.sqrt((gamma-1)*(H_hatR[m,n]-0.5*(u_hat[m,n]**2+v_hat[m,n]**2)))
      c_hatZ[m,n] = np.sqrt((gamma-1)*(H_hatZ[m,n]-0.5*(u_hat[m,n]**2+v_hat[m,n]**2)))
###----z subroutine
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

  for j in range(eqNum):
    for w in range(waveNum):
      for m in range(M+1): #spatial variable r
        for n in range(N+1): #spatial variable z
          Wg[w,j,m,n] = ag[w,m,n] * Rg[j,w,m,n]

  amdqg[:,:,:] = 0.
  apdqg[:,:,:] = 0.
 
#      print('sg = ' + str(sg))
#      print('Wg = ' + str(Wg))
   
  for j in range(eqNum): 
    for m in range(M-1):
      for n in range(N-1):
        for w in range(waveNum):
          amdqg[j,m,n] += min(sg[w,m,n+1],0) * Wg[w,j,m,n+1]
          apdqg[j,m,n] += max(sg[w,m,n],0) * Wg[w,j,m,n]   
  

  for j in range(eqNum): 
    for m in range(M+1): #q = q.shape[0],q.shape[1],q.shape[2] (eqnum,n+2,m+2)
      for n in range(N+1):
        q_new_b[j,m,n+1] = q_new_a[j,m,n+1] - (1./2.) * dt/dz * (amdqg[j,m,n] + apdqf[j,m,n])  
    
#  print('q_new_b = ' + str(q_new_b))
                ##Add ghost cell filling function
#          q_new_b = fl.Mat_ghostCells(q_new_b, M, N, 'extrap', DIM)
  print('q_new_a - q_new_b = ' + str(q_new_a - q_new_b))
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
        


###___PLOTTING___###

  if(MOVIE == 0):  
#    if(t%75 == 0 or t == 1): 
    if(t%1 == 0 or t == 1): 

      plt.figure(1)
      R,Z = np.meshgrid(rc, zc)
      plt.title(' time  = ' + str(dt*t) )
#      plt.contourf(R, Z, 
#      print('Pressure = ' + str(P))
#      print('rho = ' + str(rho))
#      print('u = ' + str(u))
#      print('v = ' + str(v))
      cp = plt.contourf(rho, label = 'density')
      plt.colorbar(cp)
#      plt.plot(rc, q_new[0, :len(rc)], label = 'hnew')
#      plt.plot(rc, q_new[1, :len(rc)], label = 'hunew') # / q[0, :len(rc)])
#      plt.plot(rc, q[0, :len(rc)], label = 'density')
#      plt.plot(rc, u[:len(rc)], label = 'velocity')# / q[0, :len(rc)])
#      plt.plot(zc, u[:len(rc)], label = 'velocity')# / q[0, :len(rc)])
#      plt.plot(rc, P[:len(rc)], label = 'pressure')# / q[0, :len(rc)])
      plt.axvline(x=0., color='k', linestyle='--')  #(max(rc)/2.))
      plt.legend()

      plt.figure(1)
#      plt.plot(rc, Q_heat[:len(rc)])
#      plt.plot(rc, EfOld[:len(rc)], label = 'EF')
#      plt.plot(rc, JOld[:len(rc)], label = 'Current Density')
#      plt.plot(rc, Temp[1:], label = 'Temperature')
      plt.legend()

#      plt.figure(3)
#      plt.plot(rc, Ndens[1:], label = 'NeutralDensity')
#      plt.legend()

      plt.show()

  if(MOVIE == 1):  
    plt.figure(1)
    plt.clf()

    plt.title(' Step number  = ' + str(t) )

    if (HEAT == 1):
      plt.plot(rc, Q_heat[:len(rc)])
    plt.plot(rc, q[0, :len(rc)], label = 'density')
    plt.plot(rc, u[:len(rc)], label = 'velocity')# / q[0, :len(rc)])
    plt.plot(rc, P[:len(rc)], label = 'pressure')# / q[0, :len(rc)])

#    plt.axvline(x=0.)  #(max(rc)/2.))
    if(TITLE == 'constant'):
      plt.axis([rc.min(),rc.max(),-1.1,2.1])
    plt.pause(1.5)
    plt.legend()
#  print('q = ' + str(q))
#  print('----')
#  EfOld = EfNew
#  JOld = JNew
  if (HEAT == 0):
    if (CYL == 0):
      q = q_new_b
    elif (CYL == 1):
      q = q_new_c

#  print('q = ' + str(q))
#  print('----')







