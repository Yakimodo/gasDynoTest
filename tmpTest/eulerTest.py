###====This is a test euler program, if broken try eulerEqs.py

import numpy as np
import matplotlib.pyplot as plt
import function_list_Test as fl

EFIX = 0
MOVIE = 0
HEAT = 0

N = 7
L = .005 #domain_length
CFL = 0.6 #Courant-Fredrichs-Lewy condition
dx = L / (2.*(N-1)) #spatial resolution
dt = CFL*dx #time step
T_begin = 1 #step number [important in determining xi = x/t ----> (in this case) --> xc/(t*dt)]
T_end = 500 #number of steps
gamma = 1.4 #gravity
xc = np.linspace(-L/2, L/2, N)

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


q, u, P = fl.initialConditions('gasDyno', 'shockTube', xc, N, q_zeros)
#q, u, P = fl.initialConditions('gasDyno', 'twoShock', xc, N, q_zeros)
#q, u, P = fl.initialConditions('gasDyno', 'flat', xc, N, q_zeros)
q_new_a = np.zeros((q.shape[0],q.shape[1]))
q_new_b = np.zeros((q.shape[0],q.shape[1]))
#  q0 = density
#  q1 = momentum
#  q2 = energy density

for t in range(T_begin, T_end+1):

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
   
    print('a = \n' + str(a))
    print('R = \n' + str(R))
    print('W = \n' + str(W))
    
    amdq[:,:] = 0.
    apdq[:,:] = 0.
 
    for j in range(eqNum): 
      for k in range(N-1):
        for w in range(waveNum):
          amdq[j,k] += min(s[w,k+1],0) * W[w,j,k+1]
          apdq[j,k] += max(s[w,k],0) * W[w,j,k]   

    print('amdq = \n' + str(amdq))
    print('apdq = \n' + str(apdq))

    if (HEAT == 0):
      for j in range(eqNum): 
        for k in range(N): #q = q.shape[0],q.shape[1] (eqnum,n+2)
          q_new_a[j,k+1] = q[j,k+1] - dt/dx * (amdq[j,k] + apdq[j,k])
       #qn+1 = qn - dt/dx*(leftFluctuations + rightFluctuations)     
      print('q_new_a = ' + str(q_new_a))

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
    if(t%75 == 0 or t == 1): 
#    if(t%5 == 0 or t == 1): 
      print('pl' + str(pl))
      print('rhol' + str(rhol))
      print('ul' + str(ul))
      print('ur' + str(ur))
      plt.figure(1)
      plt.title(' time  = ' + str(dt*t) )
#      plt.plot(xc, q_new[0, :len(xc)], label = 'hnew')
#      plt.plot(xc, q_new[1, :len(xc)], label = 'hunew') # / q[0, :len(xc)])
      plt.plot(xc, q[0, :len(xc)], label = 'density')
      plt.plot(xc, u[:len(xc)], label = 'velocity')# / q[0, :len(xc)])
      plt.plot(xc, P[:len(xc)], label = 'pressure')# / q[0, :len(xc)])
      plt.axvline(x=0., color='k', linestyle='--')  #(max(xc)/2.))
      plt.legend()
      if (HEAT == 1):
        plt.figure(2)
        plt.plot(xc, Q_heat[:len(xc)])
      plt.show()

  if(MOVIE == 1):  
    plt.figure(1)
    plt.clf()

    plt.title(' Step number  = ' + str(t) )
    plt.plot(xc, q[0, :len(xc)], label = 'density')
    plt.plot(xc, u[:len(xc)], label = 'velocity')# / q[0, :len(xc)])
    plt.plot(xc, P[:len(xc)], label = 'pressure')# / q[0, :len(xc)])
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


















