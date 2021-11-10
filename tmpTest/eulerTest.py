##====This is a test euler program, if broken try eulerEqs.py

import numpy as np
import matplotlib.pyplot as plt
import function_list_Test as fl
from matplotlib import gridspec

#Harten-Hyman Entropy Fix
EFIX = 0
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
dx = L / (N-1) #spatial resolution
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

#xo = 0.4 #position of initial discontinuity
q, u, P = fl.initialConditions('gasDyno', 'shockTube', xc, 0.3, N, q_zeros) ##ToroPg151_Test1
#q, u, P = fl.initialConditions('gasDyno', 'ToroPg151_Test2', xc, N, q_zeros)
#q, u, P = fl.initialConditions('gasDyno', 'ToroPg151_Test3', xc, 0.5, N, q_zeros)
#q, u, P = fl.initialConditions('gasDyno', 'ToroPg151_Test4', xc, 0.4, N, q_zeros)
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
  if (t%20 == 0):
    print('\nfor stepnumber = ' + str(t)) 

  q = fl.Mat_ghostCells(q, N, 'extrap' )
#  q_new_a = fl.Mat_ghostCells(q_new_a, N, 'extrap' )
#  q_new_b = fl.Mat_ghostCells(q_new_b, N, 'extrap' )
  
  if (HEAT == 1):
    Q_heat = []
    for i in range(N):
      Q_heat.append( 0.000001*np.exp(0.05*t)*np.exp(-np.abs((xc[i]-xc[np.int((N-1)/2)]))**2 / 0.002) )

  dq1, dq2, dq3 = fl.JumpSplit(q, N, 'gasDyno') ##length of dq* = N + 1
  rhol = q[0,:-1] #length of l's and r's == N+1
  rhor = q[0,1:]
  rho = q[0,:]
  moml = q[1,:-1] #rhol * ul
  momr = q[1,1:] #rhor * ur
  ul = moml/rhol #u[1,:-1]
  ur = momr/rhor#u[1,1:]
  u = q[1,:]/q[0,:] ##length == 103
  El = q[2,:-1]#pl/(gamma - 1) + 0.5 * rhol * ul**2 #q[1,:-1] / q[0,:-1]
  Er = q[2,1:]#pr/(gamma - 1) + 0.5 * rhor * ur**2 #q[1,1:] / q[0,1:]  
  pl = (gamma-1)*(El-0.5*rhol*ul**2)#P[2,:-1]
  pr = (gamma-1)*(Er-0.5*rhor*ur**2)#P[2,:-1]
  P = (gamma-1)*(q[2,:]-0.5*q[0,:]*u**2) ##length = 103
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
#    if (k == N):
#      print('a[0,:] = ' + str(a[0,:]))
  
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

  
###___Harten-Hyman Entropy Fix___###
  if (EFIX == 1):
    ll = np.zeros((2, N+1))
 #   sl = np.zeros((2, N+1))
    lr = np.zeros((2, N+1))
 #   sr = np.zeros((2, N+1))
 #   ll = np.zeros((1,N+1))
 #   lr = np.zeros((1,N+1))
    beta = 0.#np.zeros((eqNum, N+1))
    sll = np.zeros((1,N+1))
    slr = np.zeros((1,N+1))
  
    rhoml = rhol[:] + a[0,:]
    rhomr = rhor[:] - a[2,:]
    uml = ( rhol[:]*ul[:] + a[0,:]*(u_hat[:] - c_hat[:] )) / (rhol[:] + a[0,:])
    umr = ( rhor[:]*ur[:] - a[2,:]*(u_hat[:] + c_hat[:] )) / (rhor[:] - a[2,:])
    Pml = (gamma-1)*(El[:] + a[0,:]*(H_hat[:]-u_hat[:]*c_hat[:]) - 0.5*rhoml[:]*uml[:]**2)
    Pmr = (gamma-1)*(Er[:] - a[2,:]*(H_hat[:] + u_hat[:]*c_hat[:]) - 0.5*rhomr[:]*umr[:]**2)
    cml = np.sqrt(gamma*Pml/rhoml)
    cmr = np.sqrt(gamma*Pmr/rhomr)

#    for k in range(N+1):  
    for k in range(N+1):
      ll[0,k] = u[k] - c[k]
      ll[1,k] = uml[k] - cml[k]
  
      lr[1,k] = u[k] + c[k]
      lr[0,k] = umr[k] + cmr[k]
  
    for k in range(N):
      if (ll[0,k] < 0. < ll[1,k]):
        transonic = 1
        print('\nleft transonic')
        print('lkl = ' + str(ll[0,k]))
        print('lkr = ' + str(ll[1,k]))
  
      if (lr[0,k] < 0 < lr[0,k+1]):
        transonic = 2
        print('\nright transonic')
  
  
  #    transonic = 0
    tran_pos_row = 0
    tran_pos_col = 0
    if(transonic == 1):
      for k in range(N):
        if(ll[0,k] < 0. < ll[1,k]):
  #           transonic += 1
          beta = (ll[1,k] - s[0,k]) / (ll[1,k] - ll[0,k])
  #          print('\ns[0,k and k+1] before HH = ' + str(s[0,k]) + '\t' +str(s[0,k+1]))
          sll[0,:] = beta*ll[0,:]
          slr[0,:] = (1-beta)*ll[1,:]
          print('sll = ' + str(sll))
          print('slr = ' + str(slr))
  
          tran_pos_row = 0
  #        tran_pos_col = k



#    sl1, sl2 = HH_EFIX()
#    print('YES HH EFIX')
###______________________________###

  for j in range(eqNum):
    for w in range(waveNum):
      for k in range(N+1): #spatial variable
        W[w,j,k] = a[w,k] * R[w,j,k]

  amdq[:,:] = 0.
  apdq[:,:] = 0.

#  sMin = np.zeros((eqNum, N+1)) #same dims as s (speeds)
#  sMax = np.zeros((eqNum, N+1)) #same dims as s (speeds)
#  for w in range(waveNum):
#    sMin[w,:] = min(s[w,:].min(),0.)
#    sMax[w,:] = max(s[w,:].max(),0.)

  for j in range(eqNum): 
    for k in range(N):
      for w in range(waveNum):
        if (EFIX == 1 and transonic == 1 and w == tran_pos_row):# and k == tran_pos_col):
#          amdq[j,k] += sl[1,k] * W[w,j,k+1] #k+1
#          apdq[j,k] += sl[0,k] * W[w,j,k] #k
          amdq[j,k] += sll[w,k] * W[w,j,k+1] #k+1
          apdq[j,k] += slr[w,k] * W[w,j,k] #k+1
#          apdq[j,k] += sMax[w,k] * W[w,j,k] #k
        else:
#          amdq[j,k] += sMin[w,k] * W[w,j,k+1]
#          apdq[j,k] += sMax[w,k] * W[w,j,k]  
          amdq[j,k] += min(s[w,k+1],0.) * W[w,j,k+1]
          apdq[j,k] += max(s[w,k],0.) * W[w,j,k] 
#          if (s[w,k] < 0.):
#            amdq[j,k] += s[w,k] * W[w,j,k+1]
#          else:
#            apdq[j,k] += s[w,k] * W[w,j,k]


  S_abs = np.absolute(s)
#  CFL = 1. #Courant-Fredrichs-Lewy condition
  if (abs(min(s.min(),0)) > max(s.max(),0)):
    dt = CFL*dx/max(S_abs.max(),0.) #time step
  else:
    dt = CFL*dx/max(s.max(),0.)

  if (t%100 == 0):
    print('\ntime = ' + str(t*dt) + '\n step number = ' + str(t) + '\tdt = ' + str(dt))
#  print('\nstep number = ' +str(t) + '  max speed = ' + str(dx/dt))
#  print('\ntime = ' + str(t*dt) + '\n step number = ' + str(t) + '\tdt = ' + str(dt))
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
#    q_new_a = fl.Mat_ghostCells(q_new_a, N, 'extrap' )
    q_new_b = q_new_a
    for k in range(N): #q = q.shape[0],q.shape[1] (eqNum,N+2)
      q_new_b[2,k+1] = q_new_a[2,k+1] + 0.5 * dt * Q_heat[k] 

###___PLOTTING___###

  if(MOVIE == 0):  
#    if(t%75 == 0 or t == 1): 
#    if(t%5 == 0 or t == 1): 
#    if(0.249 < t*dt < 0.251 or t == 1): ##FOR Toro_Test1---Shock Tube
    if(0.195 < t*dt < 0.205 or t == 1):# or t%100 == 0): ##FOR Toro_Test1 GodunovMethod---Shock Tube 
#    if(0.148 < t*dt < 0.152 or t == 1 or t == 10 ): ##FOR Toro_Test2
#    if(0.011 < t*dt < 0.013 or t == 1 or t%10 == 0): ##FOR Toro_Test3
#    if(0.034 < t*dt < 0.036 or t == 1 or t%10 == 0): ##FOR Toro_Test4 or Toro_Test5

      print('\ntime = ' + str(t*dt) + '\n step number = ' + str(t) + '\tdt = ' + str(dt))
####____Subplots___###
#      fig,((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)#, figsize=(5, 5))
##      ax = plt.gca()
##      #sets the ratio to 1
##      ax.set_aspect(1)
#      ax1.scatter(xc, q[0, :len(xc)], label = 'density')
#      ax2.scatter(xc, u[:len(xc)], label = 'velocity')# / q[0, :len(xc)])
#      ax3.scatter(xc, P[:len(xc)], label = 'pressure')# / q[0, :len(xc)])
#      ax4.scatter(xc, e[:len(xc)], label = 'internal energy')
#      plt.axvline(x=0.3, color='k', linestyle='--')  #(max(xc)/2.))

      plt.suptitle(' time = '  + str(dt*t) + '\t timestep = '+ str(dt) )

      plt.subplot(2,2,1)
      ax = plt.gca()
      #sets the ratio to 1
      ax.set_aspect(1)
      plt.title('Mass Density' )
#      plt.plot(xc, q[0, :len(xc)], label = 'density', linestyle = 'dotted')
      plt.scatter(xc, q[0, :len(xc)], label = 'density')
      plt.xlim(0,1)
      plt.ylim(0,1.1)
      plt.axvline(x=0.3, color='k', linestyle='--')  #(max(xc)/2.))

      plt.subplot(2,2,2)
      ax = plt.gca()
      #sets the ratio to 1
#      ax.set_aspect(1./1.8)
      plt.title('Velocity')
#      plt.plot(xc, u[:len(xc)], label = 'velocity', linestyle = 'dotted')# / q[0, :len(xc)])
      plt.scatter(xc, u[:len(xc)], label = 'velocity')# / q[0, :len(xc)])
      plt.xlim(0,1)
      plt.ylim(-0.1,1.7)
      plt.axvline(x=0.3, color='k', linestyle='--')  #(max(xc)/2.))

      plt.subplot(2,2,3)
      ax = plt.gca()
      #sets the ratio to 1
      ax.set_aspect(1)
      plt.title('Pressure')
#      plt.plot(xc, P[:len(xc)], label = 'pressure', linestyle = 'dotted')# / q[0, :len(xc)])
      plt.scatter(xc, P[:len(xc)], label = 'pressure')# / q[0, :len(xc)])
      plt.xlim(0,1)
#      plt.ylim(0,1)
      plt.axvline(x=0.3, color='k', linestyle='--')  #(max(xc)/2.))

      plt.subplot(2,2,4)
      ax = plt.gca()
#      #sets the ratio to 1
      ax.set_aspect(1./2.)
      plt.title('Internal Energy')
#      plt.plot(xc, e[:len(xc)], label = 'internal energy', linestyle = 'dotted')
      plt.scatter(xc, e[:len(xc)], label = 'internal energy')
      plt.xlim(0,1)
      plt.ylim(1.8,3.8)
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


















