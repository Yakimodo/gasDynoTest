import numpy as np
import matplotlib.pyplot as plt
import function_list as fl

EFIX = 0
MOVIE = 0

N = 1001
L = 10. #domain_length
CFL = 0.6 #Courant-Fredrichs-Lewy condition
dx = L / (N-1) #spatial resolution
dt = CFL*dx #time step
T_begin = 1 #step number [important in determining xi = x/t ----> (in this case) --> xc/(t*dt)]
T_end = 500 #number of steps
gamma = 1.4 #gravity
xc = np.linspace(-L/2, L/2, N)

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

q = fl.initialConditions('gasDyno', 'shockTube', xc, N, q_empty)


for t in range(T_begin, T_end+1):

  q_new = np.empty((q.shape[0],q.shape[1]))
  q = fl.Mat_ghostCells(q, N, 'extrap')

  dq1, dq2, dq3 = fl.JumpSplit(q, N, 'gasDyno')
  rhol = q[0,:-1]
  rhor = q[0,1:]
  ul = q[1,:-1]
  ur = q[1,1:]
  pl = q[2,:-1]
  pr = q[2,1:]
  El = pl/(gamma - 1) + 0.5 * rhol * ul**2 #q[1,:-1] / q[0,:-1]
  Er = pr/(gamma - 1) + 0.5 * rhor * ur**2 #q[1,1:] / q[0,1:]  
  cl = np.sqrt(gamma*pl/rhol)
  cr = np.sqrt(gamma*pr/rhor)

  for k in range(N+1):
    u_hat[k] = (np.sqrt(rhol[k])*ul[k] + np.sqrt(rhor[k])*ur[k]) / (np.sqrt(rhol[k]) + np.sqrt(rhor[k]))
    H_hat[k] = ((El[k] + pl[k]) / np.sqrt(rhol[k]) + (Er[k] + pr[k]) / np.sqrt(rhor[k])) / (np.sqrt(rhol[k]) + np.sqrt(rhor[k])) 
    c_hat[k] = np.sqrt((gamma-1)*(H_hat[k]-0.5*u_hat[k]**2))

  
    a[1,k] = ((gamma-1)*((H_hat[k]-u_hat[k]**2)*dq1[k] + u_hat[k]*dq2[k]-dq3[k])) / c_hat[k]**2
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
    R[2,0,k] = H_hat[k] + u_hat[k]*c_hat[k]


  if (EFIX == 0):
    for j in range(eqNum):
      for w in range(waveNum):
        for k in range(N+1): #spatial variable
          W[w,j,k] = a[w,k] * R[w,j,k]
    
    amdq[:,:] = 0.
    apdq[:,:] = 0.
 
    for j in range(eqNum): 
      for k in range(0,N-1):      #____ q ==== q.shape[0],q.shape[1] ==== (eqNum,N+2) ____#
        for w in range(waveNum):
          amdq[j,k] += min(s[w,k+1],0) * W[w,j,k+1]
          apdq[j,k] += max(s[w,k],0) * W[w,j,k]   

    for j in range(eqNum): 
      for k in range(N):      #q = q.shape[0],q.shape[1] (eqNum,N+2)
        q_new[j,k+1] = q[j,k+1] - dt/dx * (amdq[j,k] + apdq[j,k])    #qn - dt/dx*(leftFluctuations + rightFluctuations)    


###___PLOTTING___###

  if(MOVIE == 0):  
    if(t%15 == 0 or t == 1):  
      plt.title(' time  = ' + str(dt*t) )
#      plt.plot(xc, q_new[0, :len(xc)], label = 'hnew')
#      plt.plot(xc, q_new[1, :len(xc)], label = 'hunew') # / q[0, :len(xc)])
      plt.plot(xc, q[0, :len(xc)], label = 'density')
      plt.plot(xc, q[1, :len(xc)], label = 'velocity')# / q[0, :len(xc)])
      plt.plot(xc, q[2, :len(xc)], label = 'pressure')# / q[0, :len(xc)])
      plt.axvline(x=0.)  #(max(xc)/2.))
      plt.legend()
      plt.show()
  if(MOVIE == 1):  
    plt.figure(1)
    plt.clf()

    plt.title(' Step number  = ' + str(t) )
    plt.plot(xc, q_new[0, :len(xc)], label = 'hnew')
    plt.plot(xc, q_new[1, :len(xc)], label = 'hunew') # / q[0, :len(xc)])
    plt.axvline(x=0.)  #(max(xc)/2.))
    plt.axis([xc.min(),xc.max(),-2,4])
    plt.pause(0.09)
#    plt.legend()
#  print('q = ' + str(q))
#  print('----')
  q = q_new
#  print('q = ' + str(q))
#  print('----')


















