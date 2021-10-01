###====This is a test euler program, if broken try eulerEqs.py

import numpy as np
import matplotlib.pyplot as plt
import function_list_Test as fl

#Harten-Hyman Entropy Fix
EFIX = 1
if (EFIX == 0):
  print('NO EFIX')
else:
  print('YES EFIX')
# To plot a movie or not
MOVIE = 0
#To add a heat source
HEAT = 0

N = 101
L = 1. #domain_length
dx = L / ((N-1)) #spatial resolution
CFL = 0.9 #Courant-Fredrichs-Lewy condition
#dt = CFL*dx #time step
T_begin = 1 #step number [important in determining xi = x/t ----> (in this case) --> xc/(t*dt)]
T_end = 8000 #number of steps
gamma = 1.4 #gravity
#xc = np.linspace(-L/2, L/2, N)
xc = np.linspace(0, L, N)

eqNum = 3 #number of equations in system (SW = 2)
waveNum = eqNum #number of waves per Riemann Prob

q_zeros = np.zeros((eqNum, N+2))
W = np.zeros((eqNum, waveNum, N+1))
R = np.zeros((eqNum, waveNum, N+1))
u_hat = np.zeros((N+1))
H_hat = np.zeros((N+1))
c_hat = np.zeros((N+1))
a = np.zeros((eqNum, N+1))
s = np.zeros((eqNum, N+1))
amdq = np.zeros((eqNum, N+1))
apdq = np.zeros((eqNum, N+1))


q, u, P = fl.initialConditions('gasDyno', 'shockTube', xc, N, q_zeros) ##ToroPg151_Test1
#q, u, P = fl.initialConditions('gasDyno', 'ToroPg151_Test2', xc, N, q_zeros)
#q, u, P = fl.initialConditions('gasDyno', 'ToroPg151_Test3', xc, N, q_zeros)
#q, u, P = fl.initialConditions('gasDyno', 'ToroPg151_Test4', xc, N, q_zeros)
#q, u, P = fl.initialConditions('gasDyno', 'ToroPg151_Test5', xc, N, q_zeros)
#q, u, P = fl.initialConditions('gasDyno', 'BalbasICs', xc, N, q_zeros)
#q, u, P = fl.initialConditions('gasDyno', 'twoShock', xc, N, q_zeros)
#q, u, P = fl.initialConditions('gasDyno', 'flat', xc, N, q_zeros)
q_new_a = np.zeros((q.shape[0],q.shape[1]))
q_new_b = np.zeros((q.shape[0],q.shape[1]))
#  q0 = density
#  q1 = momentum
#  q2 = energy density

for t in range(T_begin, T_end+1):
  if (t%100 == 0):
    print('\ntime = ' + str(t*dt) + '\n step number = ' + str(t) + '\tdt = ' + str(dt))

  q = fl.Mat_ghostCells(q, N, 'extrap' )
#  q_new_a = fl.Mat_ghostCells(q_new_a, N, 'extrap' )
#  q_new_b = fl.Mat_ghostCells(q_new_b, N, 'extrap' )
  
  if (HEAT == 1):
    Q_heat = []
    for i in range(N):
      Q_heat.append( 0.000001*np.exp(0.05*t)*np.exp(-np.abs((xc[i]-xc[np.int((N-1)/2)]))**2 / 0.002) )

  dq1, dq2, dq3 = fl.JumpSplit(q, N, 'gasDyno')
  rhol = q[0,:-1]
  rhor = q[0,1:]
  rho = q[0,:]
  moml = q[1,:-1] #rhol * ul
  momr = q[1,1:] #rhor * ur
  ul = moml/rhol #u[1,:-1]
  ur = momr/rhor#u[1,1:]
  u = q[1,:]/q[0,:]
  El = q[2,:-1]#pl/(gamma - 1) + 0.5 * rhol * ul**2 #q[1,:-1] / q[0,:-1]
  Er = q[2,1:]#pr/(gamma - 1) + 0.5 * rhor * ur**2 #q[1,1:] / q[0,1:]  
  pl = (gamma-1)*(El-0.5*rhol*ul**2)#P[2,:-1]
  pr = (gamma-1)*(Er-0.5*rhor*ur**2)#P[2,:-1]
  P = (gamma-1)*(q[2,:]-0.5*q[0,:]*u**2)
  cl = np.sqrt(gamma*pl/rhol) ##necessary for EFIX ---
  cr = np.sqrt(gamma*pr/rhor) ##----------------------
  c = np.sqrt(gamma*P/rho)
  e = P / ((gamma - 1)*q[0,:])
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

    S_abs = np.absolute(s)
    CFL = 1. #Courant-Fredrichs-Lewy condition
    if (abs(min(s.min(),0)) > max(s.max(),0)):
      dt = CFL*dx/max(S_abs.max(),0.) #time step
    else:
      dt = CFL*dx/max(s.max(),0.)

#    print('\nstep number = ' +str(t) + '  max speed = ' + str(dx/dt))
#    print('\ntime = ' + str(t*dt) + '\n step number = ' + str(t) + '\tdt = ' + str(dt))
    if (HEAT == 0):
      for j in range(eqNum): 
        for k in range(N): #q = q.shape[0],q.shape[1] (eqnum,n+2)
          q_new_a[j,k+1] = q[j,k+1] - dt/dx * (amdq[j,k] + apdq[j,k])
       #qn+1 = qn - dt/dx*(leftFluctuations + rightFluctuations)     
#      print('q_new_a = ' + str(q_new_a))

#    q_new_a = fl.Mat_ghostCells(q_new_a, N, 'extrap')

#    print(q_new_a.shape[0],q_new_a.shape[1], q_new_b.shape[0],q_new_b.shape[1],len(Q_heat))
    if (HEAT == 1):
      for j in range(eqNum): 
        for k in range(N): #q = q.shape[0],q.shape[1] (eqnum,n+2)
          q_new_a[j,k+1] = q[j,k+1] - 0.5 * dt/dx * (amdq[j,k] + apdq[j,k])
#      q_new_a = fl.Mat_ghostCells(q_new_a, N, 'extrap' )
      q_new_b = q_new_a
      for k in range(N): #q = q.shape[0],q.shape[1] (eqNum,N+2)
        q_new_b[2,k+1] = q_new_a[2,k+1] + 0.5 * dt * Q_heat[k] 


###___Harten-Hyman Entropy Fix___###
  if (EFIX == 1):
    l = np.zeros((eqNum, N+1))
    beta = np.zeros((eqNum, N+1))
    for k in range(N+1):
      l[0,k] = u[k] - c[k]
      l[1,k] = u[k]
      l[2,k] = u[k] + c[k]

    transonic = 0
    tran_pos_row = 0
    tran_pos_col = 0
    for w in range(waveNum):
      for k in range(N):
        if(l[w,k] < 0. < l[w,k+1]):
          transonic += 1
          beta[w,k] = (s[w,k] - l[w,k+1]) / (l[w,k+1] - l[w,k])
          s[w,k+1] = beta[w,k]*l[w,k]
          s[w,k] = (1-beta[w,k])*l[w,k+1]
          tran_pos_row = w
          tran_pos_col = k

    for j in range(eqNum):
      for w in range(waveNum):
        for k in range(N+1): #spatial variable
          W[w,j,k] = a[w,k] * R[w,j,k]
    
    amdq[:,:] = 0.
    apdq[:,:] = 0.   
 
    for j in range(eqNum): 
      for k in range(N-1):
        for w in range(waveNum):
          if (transonic > 0 and w == tran_pos_row and k == tran_pos_col):
            amdq[j,k] += s[w,k+1] * W[w,j,k+1]
            apdq[j,k] += s[w,k] * W[w,j,k]   
          else:
            amdq[j,k] += min(s[w,k+1],0) * W[w,j,k+1]
            apdq[j,k] += max(s[w,k],0) * W[w,j,k]   

    S_abs = np.absolute(s)
    CFL = 1. #Courant-Fredrichs-Lewy condition
    if (abs(min(s.min(),0)) > max(s.max(),0)):
      dt = CFL*dx/max(S_abs.max(),0.) #time step
    else:
      dt = CFL*dx/max(s.max(),0.)

    print('\nstep number = ' +str(t) + '  max speed = ' + str(dx/dt))
    print('\ntime = ' + str(t*dt) + '\n step number = ' + str(t) + '\tdt = ' + str(dt))
    if (HEAT == 0):
      for j in range(eqNum): 
        for k in range(N): #q = q.shape[0],q.shape[1] (eqnum,n+2)
          q_new_a[j,k+1] = q[j,k+1] - dt/dx * (amdq[j,k] + apdq[j,k])
       #qn+1 = qn - dt/dx*(leftFluctuations + rightFluctuations)     
#      print('q_new_a = ' + str(q_new_a))

#    q_new_a = fl.Mat_ghostCells(q_new_a, N, 'extrap')

#    print(q_new_a.shape[0],q_new_a.shape[1], q_new_b.shape[0],q_new_b.shape[1],len(Q_heat))
    if (HEAT == 1):
      for j in range(eqNum): 
        for k in range(N): #q = q.shape[0],q.shape[1] (eqnum,n+2)
          q_new_a[j,k+1] = q[j,k+1] - 0.5 * dt/dx * (amdq[j,k] + apdq[j,k])
#      q_new_a = fl.Mat_ghostCells(q_new_a, N, 'extrap' )
      q_new_b = q_new_a
      for k in range(N): #q = q.shape[0],q.shape[1] (eqNum,N+2)
        q_new_b[2,k+1] = q_new_a[2,k+1] + 0.5 * dt * Q_heat[k] 

###___PLOTTING___###

  if(MOVIE == 0):  
#    if(t%75 == 0 or t == 1): 
#    if(t%5 == 0 or t == 1): 
#    if(0.249 < t*dt < 0.251 or t == 1): ##FOR Toro_Test1---Shock Tube
    if(0.199 < t*dt < 0.201 or t == 1 or t%100 == 0): ##FOR Toro_Test1 GodunovMethod---Shock Tube 
#    if(0.148 < t*dt < 0.152 or t == 1 or t == 10 ): ##FOR Toro_Test2
#    if(0.01198 < t*dt < 0.0121 or t == 1): ##FOR Toro_Test3
#    if(0.0349 < t*dt < 0.0351 or t == 1): ##FOR Toro_Test4 or Toro_Test5
      print('\ntime = ' + str(t*dt) + '\n step number = ' + str(t) + '\tdt = ' + str(dt))
      plt.suptitle(' time = '  + str(dt*t) + '\t timestep = '+ str(dt) )
      plt.subplot(2,2,1)
      plt.title('Mass Density' )
      plt.plot(xc, q[0, :len(xc)], label = 'density', linestyle = 'dotted')
      plt.axvline(x=0.3, color='k', linestyle='--')  #(max(xc)/2.))

      plt.subplot(2,2,2)
      plt.title('Velocity')
      plt.plot(xc, u[:len(xc)], label = 'velocity', linestyle = 'dotted')# / q[0, :len(xc)])
      plt.axvline(x=0.3, color='k', linestyle='--')  #(max(xc)/2.))

      plt.subplot(2,2,3)
      plt.title('Pressure')
      plt.plot(xc, P[:len(xc)], label = 'pressure', linestyle = 'dotted')# / q[0, :len(xc)])
      plt.axvline(x=0.3, color='k', linestyle='--')  #(max(xc)/2.))

      plt.subplot(2,2,4)
      plt.title('Internal Energy')
      plt.plot(xc, e[:len(xc)], label = 'internal energy', linestyle = 'dotted')
      plt.axvline(x=0.3, color='k', linestyle='--')  #(max(xc)/2.))
#      plt.legend()
#    plt.title('mass density' )
#    cp = plt.contourf(rho[:,:])
#    plt.colorbar(cp)
#    plt.subplot(2,2,2)
#    plt.title('Pressure' )
#    cp = plt.contourf(P[:,:])
#    plt.colorbar(cp)
    plt.show()

#      print('pl' + str(pl))
#      print('rhol' + str(rhol))
#      print('ul' + str(ul))
#      print('ur' + str(ur))
#      plt.figure(1)
#      plt.title(' time and timestep = ' + str(dt*t) + str(dt) )
##      plt.plot(xc, q_new[0, :len(xc)], label = 'hnew')
##      plt.plot(xc, q_new[1, :len(xc)], label = 'hunew') # / q[0, :len(xc)])
#      plt.plot(xc, q[0, :len(xc)], label = 'density')
#      plt.plot(xc, u[:len(xc)], label = 'velocity')# / q[0, :len(xc)])
#      plt.plot(xc, P[:len(xc)], label = 'pressure')# / q[0, :len(xc)])
#      plt.plot(xc, e[:len(xc)], label = 'internal energy')
#      plt.axvline(x=0., color='k', linestyle='--')  #(max(xc)/2.))
#      plt.legend()
#      if (HEAT == 1):
#        plt.figure(2)
#        plt.plot(xc, Q_heat[:len(xc)])
#      plt.show()

  if(MOVIE == 1):  
    plt.figure(1)
    plt.clf()

    plt.title(' Step number  = ' + str(t) )
    plt.plot(xc, q[0, :len(xc)], label = 'density')
    plt.plot(xc, u[:len(xc)], label = 'velocity')# / q[0, :len(xc)])
    plt.plot(xc, P[:len(xc)], label = 'pressure')# / q[0, :len(xc)])
    plt.plot(xc, e[:len(xc)], label = 'internal energy')
#    plt.plot(xc, q_new[0, :len(xc)], label = 'hnew')
#    plt.plot(xc, q_new[1, :len(xc)], label = 'hunew') # / q[0, :len(xc)])
    plt.axvline(x=0.)  #(max(xc)/2.))
#    plt.axis([xc.min(),xc.max(),-2,4])
    plt.pause(0.09)
    plt.legend()
#  print('q = ' + str(q))
#  print('----')
  if (HEAT == 0):
    q = q_new_a
  elif (HEAT == 1):
    q = q_new_b
#  print('q = ' + str(q))
#  print('----')


















