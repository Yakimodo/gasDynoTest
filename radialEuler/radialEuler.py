###====This is a test euler program, if broken try eulerEqs.py

import numpy as np
import matplotlib.pyplot as plt
import function_list_Test as fl

EFIX = 0
MOVIE = 1
HEAT = 0
RADIAL = 0
TITLE = 'shockTube'
#TITLE = 'constant'
#TITLE = 'twoShock'


N = 100
L = 10. #domain_length
CFL = 0.6 #Courant-Fredrichs-Lewy condition
dr = L / (2.*(N-1)) #spatial resolution
dt = CFL*dr #time step
T_begin = 1 #step number [important in determining xi = x/t ----> (in this case) --> xc/(t*dt)]
T_end = 500 #number of steps
gamma = 1.4 #gravity
kB = 1.38e-23
rc = np.linspace(0., L, N+1)


eqNum = 3 #number of equations in system (SW = 2)
waveNum = eqNum #number of waves per Riemann Prob

q_empty = np.empty((eqNum, N+2))
W = np.empty((eqNum, waveNum, N+1))
R = np.empty((eqNum, waveNum, N+1))
u_hat = np.empty((N+1))
H_hat = np.empty((N+1))
c_hat = np.empty((N+1))
a = np.empty((eqNum, N+1))
s = np.empty((eqNum, N+1))
amdq = np.empty((eqNum, N+1))
apdq = np.empty((eqNum, N+1))
Q_heat = []

q, u, P = fl.initialConditions('gasDyno', TITLE, rc, N, q_empty)
#q, u, P = fl.initialConditions('gasDyno', 'constant', rc, N, q_empty)
#q, u, P = fl.initialConditions('gasDyno', 'twoShock', rc, N, q_empty)
#  q0 = density
#  q1 = velocity
#  q2 = pressure

Ndens = []
for j in range(N):
  Ndens.append(1.e20*np.exp(-np.abs((rc[j]-rc[np.int((N-1)/2)]))**2 / 0.2))

Ndens = fl.fill_ghosts(Ndens, N, 'extrap')
Ndens = np.array(Ndens)


for t in range(T_begin, T_end+1):

  q_new_a = np.empty((q.shape[0],q.shape[1]))
  q_new_b = np.empty((q.shape[0],q.shape[1]))

  if (RADIAL == 0):
    q = fl.Mat_ghostCells(q, N, 'extrap')
  if (RADIAL == 1):
    q = fl.Mat_ghostCells(q, N, 'wall')

  EfOld = []
  for j in range(N):
    EfOld.append(np.exp(-np.abs((rc[j]-rc[np.int((N-1)/4)]))**2 / 0.02))
  
  JOld = []
  for j in range(N):
    JOld.append(np.exp(-np.abs((rc[j]-rc[np.int((N-1)/4)]))**2 / 0.02))
  
  JOld = fl.fill_ghosts(JOld, N, 'extrap')
  EfOld = fl.fill_ghosts(EfOld, N, 'extrap') 

  Q_heat = []
  if (HEAT == 1):
    for i in range(N+1):
      Q_heat.append(25.*np.exp(-0.05*t)*EfOld[i]*JOld[i])
#      Q_heat.append( 0.01*np.exp(0.5*t)*np.exp(-np.abs((rc[i]-rc[(N-1)/2]))**2 / 0.2) )
  
  c = 1.
#  EfNew = fl.simplePropRight(EfOld, dt, dr, -c)
#  JNew = fl.simplePropLeft(JOld, dt, dr, -c)

  dq1, dq2, dq3 = fl.JumpSplit(q, N, 'gasDyno')
  rhol = q[0,:-1]
  rhor = q[0,1:]
  moml = q[1,:-1] #rhol * ul
  momr = q[1,1:] #rhor * ur
  ul = moml/rhol #u[1,:-1]
  ur = momr/rhor#u[1,1:]
  u = q[1,:]/q[0,:]
  El = q[2,:-1]#pl/(gamma - 1) + 0.5 * rhol * ul**2 #q[1,:-1] / q[0,:-1]
  Er = q[2,1:]#pr/(gamma - 1) + 0.5 * rhor * ur**2 #q[1,1:] / q[0,1:]  
  pl = (gamma-1)*(El-0.5*rhol*ul**2)#P[2,:-1]
  pr = (gamma-1)*(Er-0.5*rhor*ur**2)#P[2,:-1]
#  Temp = 2./(5.*Ndens[:]*kB) * (q[2,:] - 0.5 * q[1,:]**2/q[0,:])
  P = (gamma-1)*(q[2,:]-0.5*q[0,:]*u**2)
  cl = np.sqrt(gamma*pl/rhol)
  cr = np.sqrt(gamma*pr/rhor)

  for k in range(N+1):
    u_hat[k] = (np.sqrt(rhol[k])*ul[k] + np.sqrt(rhor[k])*ur[k]) / (np.sqrt(rhol[k]) + np.sqrt(rhor[k]))
    H_hat[k] = ((El[k] + pl[k]) / np.sqrt(rhol[k]) + (Er[k] + pr[k]) / np.sqrt(rhor[k])) / (np.sqrt(rhol[k]) + np.sqrt(rhor[k])) 
    c_hat[k] = np.sqrt((gamma-1)*(H_hat[k]-0.5*u_hat[k]**2))

  
    a[1,k] = (gamma-1)*((H_hat[k]-u_hat[k]**2)*dq1[k] + u_hat[k]*dq2[k]-dq3[k]) / c_hat[k]**2
    a[2,k] = (dq2[k] + (c_hat[k] - u_hat[k])*dq1[k] - c_hat[k]*a[1,k]) / (2*c_hat[k])
    a[0,k] = dq1[k] - a[1,k] - a[2,k]
  
    s[0,k] = u_hat[k] - c_hat[k]
    s[1,k] = u_hat[k]
    s[2,k] = u_hat[k] + c_hat[k]
  
    R[0,0,k] = 1.
    R[0,1,k] = u_hat[k] - c_hat[k]
    R[0,2,k] = H_hat[k] - u_hat[k]*c_hat[k]
    R[1,0,k] = 1.
    R[1,1,k] = u_hat[k]
    R[1,2,k] = 0.5 * u_hat[k]**2
    R[2,0,k] = 1.
    R[2,1,k] = u_hat[k] + c_hat[k]
    R[2,2,k] = H_hat[k] + u_hat[k]*c_hat[k]


  if (EFIX == 0):
    for j in range(eqNum):
      for w in range(waveNum):
        for k in range(N+1): #spatial variable
          W[w,j,k] = a[w,k] * R[w,j,k]
    
    amdq[:,:] = 0.
    apdq[:,:] = 0.
 
    for j in range(eqNum): 
      for k in range(N-1):
        for w in range(waveNum):
          amdq[j,k] += min(s[w,k+1],0) * W[w,j,k+1]
          apdq[j,k] += max(s[w,k],0) * W[w,j,k]   

    if (HEAT == 0):
      if (RADIAL == 0):
        for j in range(eqNum): 
          for k in range(N): #q = q.shape[0],q.shape[1] (eqNum,N+2)
            q_new_a[j,k+1] = q[j,k+1] - dt/dr * (amdq[j,k] + apdq[j,k])   

      if (RADIAL == 1):    
        for j in range(eqNum): 
          for k in range(N): #q = q.shape[0],q.shape[1] (eqNum,N+2)
            q_new_a[j,k+1] = q[j,k+1] - 0.5*dt/dr * (amdq[j,k] + apdq[j,k])   

            if (j == 0):
              q_new_b[j,k+1] = q_new_a[j,k+1] - 0.5*dt * (1/rc[k+1] * q_new_a[1,k+1])
            elif (j == 1):
              q_new_b[j,k+1] = q_new_a[j,k+1] - 0.5*dt * (1/rc[k+1] * (q_new_a[1,k+1]**2/q_new_a[0,k+1]))
            elif (j == 2):
              q_new_b[j,k+1] = q_new_a[j,k+1] - 0.5*dt * (1/rc[k+1] * ((q_new_a[2,k+1]+P[k+1])*(q_new_a[1,k+1]/q_new_a[0,k+1])))
       
       #qn+1 = qn - dt/dr*(leftFluctuations + rightFluctuations)    

#    q_new_a = fl.Mat_ghostCells(q_new_a, N, 'extrap')

#    print(q_new_a.shape[0],q_new_a.shape[1], q_new_b.shape[0],q_new_b.shape[1],len(Q_heat))
    if (HEAT == 1):
      for j in range(eqNum): 
        for k in range(N): #q = q.shape[0],q.shape[1] (eqNum,N+2)
          q_new_a[j,k+1] = q[j,k+1] - 0.5*dt/dr * (amdq[j,k] + apdq[j,k])
      q_new_b = q_new_a

      for k in range(N-1): #q = q.shape[0],q.shape[1] (eqNum,N+2)
        q_new_b[2,k+1] = q_new_a[j,k+1] + 0.5 * dt * Q_heat[k]
###___PLOTTING___###

  if(MOVIE == 0):  
#    if(t%75 == 0 or t == 1): 
    if(t%75 == 0 or t == 1): 
      print('pl' + str(pl))
      print('rhol' + str(rhol))
      print('cl' + str(cl))
      print('cr' + str(cr))
      print('ul' + str(ul))
      print('ur' + str(ur))
      plt.figure(1)
      plt.title(' time  = ' + str(dt*t) )
#      plt.plot(rc, q_new[0, :len(rc)], label = 'hnew')
#      plt.plot(rc, q_new[1, :len(rc)], label = 'hunew') # / q[0, :len(rc)])
      plt.plot(rc, q[0, :len(rc)], label = 'density')
      plt.plot(rc, u[:len(rc)], label = 'velocity')# / q[0, :len(rc)])
      plt.plot(rc, P[:len(rc)], label = 'pressure')# / q[0, :len(rc)])
      plt.axvline(x=0., color='k', linestyle='--')  #(max(rc)/2.))
      plt.legend()

      plt.figure(1)
      plt.plot(rc, Q_heat[:len(rc)])
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
    q = q_new_a
  elif (HEAT == 1):
    q = q_new_b
#  print('q = ' + str(q))
#  print('----')







