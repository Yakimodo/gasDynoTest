import numpy as np
import matplotlib.pyplot as plt
import function_list as fl

N = 1250
L = 2.5 #domain_length
CFL = 0.9 #Courant-Fredrichs-Lewy condition
dx = L / (2.*(N-1)) #spatial resolution
dt = CFL*dx #time step
T_begin = 1 #step number [important in determining xi = x/t ----> (in this case) --> xc/(t*dt)]
T_end = 1800 #number of steps
g = 1. #gravity
#xc = np.linspace(-L/2, L/2, N)
xc = np.linspace(0., L, N+1)

eqNum = 2 #number of equations in system (SW = 2)
waveNum = eqNum #number of waves per Riemann Prob

q_empty = np.empty((eqNum, N+2))
W = np.empty((eqNum, waveNum, N+1))
R = np.empty((eqNum, waveNum, N+1))
h_avg = np.empty((N+1))
u_hat = np.empty((N+1))
c_hat = np.empty((N+1))
a = np.empty((eqNum, N+1))
s = np.empty((eqNum, N+1))
amdq = np.empty((eqNum, N+1))
apdq = np.empty((eqNum, N+1))

RADIAL = 1
MOVIE = 0
EFIX = 0

#q = fl.initialConditions('ShallowWater', 'SmallDisturbance', xc, N, q_empty)
#q = fl.initialConditions('ShallowWater', 'DamBreak', xc, N, q_empty)
q = fl.initialConditions('ShallowWater', 'RadialDamBreak', xc, N, q_empty)

for t in range(T_begin, T_end+1):

  q_new_a = np.empty((q.shape[0],q.shape[1]))
  q_new_b = np.empty((q.shape[0],q.shape[1]))
#  q = fl.Mat_ghostCells(q, N, 'extrap')
  q = fl.Mat_ghostCells(q, N, 'wall')

  dq1, dq2, dq3 = fl.JumpSplit(q, N, 'ShallowWater')
  hl = q[0,:-1]
  hr = q[0,1:]
  hul = q[1,:-1]
  hur = q[1,1:]
  ul = q[1,:-1] / q[0,:-1]
  ur = q[1,1:] / q[0,1:]  
  cl = np.sqrt(g * hl)
  cr = np.sqrt(g * hr)
 
  for k in range(N+1):
    h_avg[k] = 0.5 * (hl[k] + hr[k])
    u_hat[k] = (np.sqrt(hr[k])*ur[k] + np.sqrt(hl[k])*ul[k]) / (np.sqrt(hl[k]) + np.sqrt(hr[k]))
    c_hat[k] = np.sqrt(g * h_avg[k])
  
    a[0,k] = ((u_hat[k] + c_hat[k])*dq1[k] - dq2[k]) / (2*c_hat[k])
    a[1,k] = (-(u_hat[k] - c_hat[k])*dq1[k] + dq2[k]) / (2*c_hat[k])
  
    s[0,k] = u_hat[k] - c_hat[k]
    s[1,k] = u_hat[k] + c_hat[k]
  
    R[0,0,k] = 1.
    R[0,1,k] = u_hat[k] - c_hat[k]
    R[1,0,k] = 1.
    R[1,1,k] = u_hat[k] + c_hat[k]
  
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
        q_new_a[j,k+1] = q[j,k+1] - 0.5*dt/dx * (amdq[j,k] + apdq[j,k])    #qn - dt/dx*(leftFluctuations + rightFluctuations)    
   
    if (RADIAL == 1):
      for j in range(eqNum):
        for k in range(N):
          if (j == 0):
            q_new_b[j,k+1] = q_new_a[j,k+1] - 0.5 * dt * (1/xc[k+1] * (q_new_a[1,k+1]))
          if (j == 1):
            q_new_b[j,k+1] = q_new_a[j,k+1] - 0.5 * dt * (1/xc[k+1] * (q_new_a[1,k+1]**2/q_new_a[0,k+1]))


####_____________#####

#  if (EFIX == 1):
#    s_MAX = []
#    s_MIN = []
#    hm = hr - a2
#    hum = hur - a2 * (u_hat + c_hat)
##    hm = q[0,:-1] + a1
##    hum = q[1,:-1] + a1 * (u_hat - c_hat)
#    um = hum/hm ## m --- intermediate state
#    qm = np.array([hm,hum])
#    lambda1_m = um - np.sqrt(g*hm)
#    lambda2_m = um + np.sqrt(g*hm)
#    lambda1_l = ul - cl
#    lambda2_r = ur + cr
#    print(lambda1_m < 0 < lambda2_r)
#    print(lambda1_l < 0 < lambda1_m)
#
###assuming only 1-transonic wave
#    for i in range(N):
#      if (lambda1_l[i] < 0. < lambda1_m[i]):
#        transonic = 1
#      elif (lambda2_m[i] < 0. < lambda2_r[i]):
#        transonic = 2
#      else:
#        transonic = 0
#
#      if (transonic == 0):
#        print('HERE0')
#        S = np.min(s1)
#        S0 = np.min(s2)
#        S1 = min(S,S0)
#        s_MIN.append(min(S1,0))
#   
#        S = np.max(s1)
#        S0 = np.max(s2)
#        S1 = max(S,S0)
#        s_MAX.append(max(S1,0))
#        
#      elif (transonic == 1):
#        print('HERE1')
#        beta = (lambda1_m[i] - s1[i]) / (lambda1_m[i] - lambda1_l[i])
#        s_MIN.append(beta*lambda1_l[i]) #lambda1_hat_minus
#        s_MAX.append((1-beta)*lambda1_m[i]) #lambda1_hat_plus
#
##        a1l = beta * a1
##        alr = (1-beta) * a1
##        W1l = a1l * 1
##        W1r = a1r * 1
#    
#      elif (transonic == 2):
#        print('HERE2')
#        beta = (lambda2_r - s2) / (lambda2_r - lambda2_m)
#        s_MIN.append(beta*lambda2_m[i]) #lambda2_hat_minus
#        s_MAX.append((1-beta)*lambda2_r[i]) #lambda2_hat_plus 
#
##        a2l = beta * a2  
##        a2r = (1-beta) * a2
#  
#  
#    W[0,0,:] = a1 * 1
#    W[0,1,:] = a1 * (u_hat - c_hat)
#    W[1,0,:] = a2 * 1
#    W[1,1,:] = a2 * (u_hat + c_hat)
#
#    for j in range(eqNum): 
#      for k in range(0,q.shape[1]-2):      #q = q.shape[0],q.shape[1] (eqNum,N+2)
#        q_new[j,k+1] = q[j,k+1] - dt/dx * (s_MAX[k]*W[0,j,k] + s_MIN[k]*W[1,j,k+1])    #qn - dt/dx*(fluctuations)    
#
#  
#  xi = (xc)/(t*dt)

  if(MOVIE == 0):  
    if(t%1700 == 0 or t == 1):# %75 == 0 or t == 1):  
      print('time step = ' + str(t*dt))
      print('step number = ' + str(t))
      plt.title(' time  = ' + str(dt*t) )
#      plt.plot(xc, q_new[0, :len(xc)], label = 'hnew')
#      plt.plot(xc, q_new[1, :len(xc)], label = 'hunew') # / q[0, :len(xc)])
      plt.plot(xc, q[0, :len(xc)], label = 'h')
      plt.plot(xc, q[1, :len(xc)], label = 'hu')# / q[0, :len(xc)])
      plt.axvline(x=0.)  #(max(xc)/2.))
      plt.legend()
      plt.show()
  if(MOVIE == 1):  
    plt.figure(1)
    plt.clf()

    plt.title(' Step number  = ' + str(t) )
    plt.plot(xc, q[0, :len(xc)], label = 'h')
    plt.plot(xc, q[1, :len(xc)], label = 'hu') # / q[0, :len(xc)])
    plt.axvline(x=0.)  #(max(xc)/2.))
    plt.axis([xc.min(),xc.max(),-0.2,2.1])
    plt.pause(0.09)
#    plt.legend()
#  print('q = ' + str(q))
#  print('----')
  if (RADIAL == 1):
    q = q_new_b
  elif (RADIAL == 0):
    q = q_new_a
#  print('q = ' + str(q))
#  print('----')




