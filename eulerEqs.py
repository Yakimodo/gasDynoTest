import matplotlib.pyplot as plt
import numpy as np
import csv
import function_list as fl
from numpy import linalg as LA

#What is needed? 
# Solution to mass continuity equation --- 1D Cartesian : dp/dt = - div ( rho * vel ) 
# ICs = Density and Velocity
MOVIE = 1 # 1=yes 0=no(saves as npy)
BC = 'Continuous' #'Periodic'#'Continous' # 'Periodic'
#sim_type = 'Sound'
#sim_details = 'VImp'
sim_type = 'ShallowWater'
sim_details = 'SmallDisturbance'
#sim_details = 'DamBreak'


##Lists and arrays----
#rho_update = []
#rho_new = []
rho = []
rho_original = []
vel = []
mom = []
gradP = []
engyDens = []
temp = []
s = []
c_sound = []

##Constants----
N = 501#number of grid cells
length = 20
#ng = 2 #number of ghost cells
CFL = 0.8 #0.0001
dx = 0.01*N
dt = CFL*dx
beta = 200.
T_end = 210
xc = np.linspace(0, length, N)
kB = 1.38e-23
#Ndens = 1e21

Z1 = 0.
Z2 = 0.
a1 = 0.
a2 = 0.
W1 = np.empty((N+2, 2))
W2 = np.empty((N+2, 2))
amdq = np.empty((N+2, 2))
apdq = np.empty((N+2, 2))
##----Evolve in time----
for t in range(0,T_end):
  print(t)

##_____Initialize_____#
  if (t==0):
    u, P, q, Z, K, rho = fl.initialConditions(sim_type, sim_details, xc, N)
    Z = fl.fill_ghosts(Z, N, 'extrap')  #Sound
###_________________________________###

###_________GHOST CELLS_____________###
  q_new = np.empty((N+2,2))
  q = fl.Mat_ghostCells(q, N, 'extrap')   
  P = fl.fill_ghosts(P, N, 'extrap')
  u = fl.fill_ghosts(u, N, 'extrap')  
  c_water = np.sqrt(P)
  if (sim_type == 'ShallowWater'):
    s = c_water
  elif (sim_type == 'Sound'):
    for i in range(len(rho)):
      c_sound.append(np.sqrt(K[i]/rho[i]))
    c_sound = fl.fill_ghosts(c_sound, N, 'extrap')
    s = c_sound
###_________________________________###

###____________PLOT_________________###
  if (MOVIE == 0):
    arr1 = np.array(P) 
    arr2 = np.array(u)
    if (t == 0): 
      shape = (T_end,len(arr1))
      matrix1 = np.empty(shape)
      matrix2 = np.empty(shape)
      MAT1 = matrix1
    matrix1[t,:] = arr1
    matrix2[t,:] = arr2
    np.save("Qevolv1" + sim_type + sim_details + ".npy", matrix1)
    np.save("Qevolv2" + sim_type + sim_details + ".npy", matrix2)

  if (MOVIE == 1):
    plt.figure(1)
    plt.clf()
    plt.title(' Step number  = ' + str(t) )
    plt.plot(xc, P[1:len(xc)+1])
    plt.plot(xc, u[1:len(xc)+1])
    plt.axvline(x=0.5*max(xc))
#    plt.axis([0, max(xc), -1., 2])
    plt.pause(0.1)
    plt.legend()
#    plt.show()

    if (t == T_end-1):
      plt.show()
###_________________________________###

  c_water = np.sqrt(P)
  for i in range(1,N+1):
    dq1, dq2 = fl.JumpSplit(q, N, i)
    #print(dq1,dq2)
###_____Sound Waves____###
    if (sim_type == 'Sound'):
      Zr = Z[i]
      Zl = Z[i-1]
      a1 = (-dq1 + Zr*dq2) / (Zr + Zl)
      a2 = (dq1 + Zl*dq2) / (Zr + Zl)
      W1[i,0] = -a1*Zl
      W1[i,1] = a1
      W2[i,0] = a2*Zr
      W2[i,1] = a2
###____________________###

###_____Shallow Water Waves_____###
    if (sim_type == 'ShallowWater'):
      cr = c_water[i]
      cl = c_water[i-1]
      a1 = (cr*dq1 - dq2) / (cr+cl)
      a2 = (cl*dq1 + dq2) / (cr+cl)
      W1[i,0] = a1
      W1[i,1] = -cl*a1
      W2[i,0] = a2
      W2[i,1] = cr*a2

  for j in range(2):
    for i in range(1, N+1):
      amdq[i,j] = -s[i+1]*W1[i+1,j]
      apdq[i,j] = s[i]*W2[i,j]

 
  for j in range(2):
    for i in range(1,N+1):
      q_new[i,j] = q[i,j] - (dt/dx) * (amdq[i,j] + apdq[i,j])

  q = q_new   

  P = []
  u = []
  if (sim_type == 'Sound'):
    for i in range(1,N+1):
      P.append(q[i,0])
      u.append(q[i,1])
  if (sim_type == 'ShallowWater'):
    for i in range(1,N+1):
      P.append(q[i,0])
      u.append(q[i,1])
  
    
#plt.figure(1)
#  plt.plot(xc, P)
#  plt.show()

###---Initialize-----
#  if (t == 0):
#    rho, vel, engyDens, Pressure, rho_original = fl.initialConditions(x,N)
#    rhoU = rho
##    rhoF = rho
#    K = rho
#    Zo = 1
#    
#   R = np.array([[-Zo, Zo],[1, 1]])
#   Rminus = LA.inv(R)
#
#
###--Interface states----
#  Qr_cont = []
#  Ql_cont = []
#  Qr_mom = []
#  Ql_mom = []
#  Qr_engyDens = []
#  Ql_engyDens = []
#
####---Change velocity----###
##  if (t>0):
##    vel_new = []
##    vel_new = fl.VelChange(vel, N, t)
##    vel = vel_new
##  Pressure = fl.fill_ghosts(Pressure, ng, N, BC)
#  rho = fl.fill_ghosts(rho, ng, N, BC)
##  print(len(Pressure), len(rho))
#  if (t != 0):
#    rhoU = fl.fill_ghosts(rhoU, ng, N, BC)
##    rhoF = fl.fill_ghosts(rhoF, ng, N, BC)
#  vel = fl.fill_ghosts(vel, ng, N, BC)
##  engyDens = fl.fill_ghosts(engyDens, ng, N,BC)
#
####---initializing/advecting temperature-------###
##  if (t == 0):
##    temp = fl.Temp(Ndens, engyDens, rho, vel, N) #initialize temperature
##  temp = fl.fill_ghosts(temp, ng, N)
#
####----intializaing/advecting momentum------###
##  mom = []
##  for i in range(len(vel)):
##    mom.append(rho[i] * vel[i])
#
##  print(len(Pressure))  
##  gradP = fl.CentDiff(Pressure, dx, N, BC)
#  Qr_cont, Ql_cont = fl.states(rho, dx, dt, vel, 'SuperBee', N, BC)
##  Qr_cont, Ql_cont = fl.states(rho, dx, dt, vel, 'MC Limiter', N)
##  Qr_mom, Ql_mom = fl.states(mom, dx, dt, vel, 'SuperBee', N, BC)
##  Qr_engyDens, Ql_engyDens = fl.states(engyDens, dx, dt, vel, 'SuperBee', N, BC)
##
#
###--Solve Jump Discontinuity----
#  Flux_cont = []
##  Flux_mom = []
##  Flux_engyDens = []
#  Flux_cont = fl.Riemann(Ql_cont, Qr_cont, vel, N)
##  Flux_mom = fl.Riemann(Ql_mom, Qr_mom, vel, N)
##  Flux_engyDens = fl.Riemann(Ql_engyDens, Qr_engyDens, vel, N)
##
#  Flux_cont = fl.fill_ghosts(Flux_cont, ng, N, BC)
##  Flux_mom = fl.fill_ghosts(Flux_mom, ng, N, BC)
##  Flux_engyDens = fl.fill_ghosts(Flux_mom, ng, N, BC) 
#
###--To update solution----
#  rho_new = []
#  rho_newUpwind = []
#  rho_newFromms = []
#  mom_new = [] 
#  engyDens_new = []
##  Qt = []
##  Qt = fl.PowerDeposition(N, t)
##  Qt = fl.fill_ghosts(Qt, ng, N)
#
## print(len(mom), len(Flux_mom), len(gradP)) 
#  for i in range(1,N+1):
####_______FLUX LIMITER____#
#    rho_new.append(rho[i] + (dt/dx) * (Flux_cont[i] - Flux_cont[i+1]))
##    mom_new.append(mom[i] + (dt/dx) * (Flux_mom[i] - Flux_mom[i+1]) + gradP[i]*(dt))
##_________________________#
#
#####______UPWIND_______####
##    rho_newUpwind.append(rhoU[i] * (1 - 2*(dt/dx)*vel[i]) + (dt/dx)*(vel[i]*rhoU[i-1] + rhoU[i]*vel[i])) 
#    rho_newUpwind.append(rhoU[i] - (dt/dx) * ( vel[i]*rhoU[i] - vel[i-1]*rhoU[i-1]))
##_________________________#
##    rho_new.append(rho[i] - (dt/(2*dx))*vel[i]*(rho[i+1] - rho[i-1]))
##    rho_new.append(rho[i] - (dt/(2*dx))*vel[i]*(rho[i+1] - rho[i-1]) +  0.5*(dt/dx)**2*vel[i]**2*(rho[i-1] - 2*rho[i] + rho[i+1]))
####_____FROMM'S METHOD______###
##    rho_newFromms.append(rhoF[i] - (dt/(2*dx))*vel[i]*(3*rhoF[i] - 4*rhoF[i-1] + rhoF[i-2]) +  0.5*(dt/dx)**2*vel[i]**2*(rhoF[i] - 2*rhoF[i-1] + rhoF[i-2]))
##_____________________________#
##    engyDens_new.append(engyDens[i] - (dt/dx)*(Flux_engyDens[i] - Flux_engyDens[i+1]))# + dt*Qt[i])
# 
###______UPDATE_____##
#  rho = rho_new
#  rhoU = rho_newUpwind
##  rhoF = rho_newFromms
##  engyDens = engyDens_new
##  mom = mom_new
###_________________##
#
###---updating temperature-----##
##  temp = fl.Temp(Ndens, engyDens, rho, vel, N-1)
##
##  temp.insert(0, temp[-1])
#
###---updating pressure------##
##  Pressure_new = []
##  Pressure = []
##  for i in range(N):
###    Pressure_new.append(Ndens * kB * temp[i])
##    Pressure.append(rho[i])
##  Pressure_new.insert(0, Pressure[-1])
##  Pressure = Pressure_new
#
##  print(vel)
###---velocity matters----##
##  vel = []
##  for i in range(len(rho)):
##    vel.append(mom_new[i] / rho[i])
#  vel = vel[1:-1]
####-----PLOT-----#####
##  fl.Porta_plotty2(x, rho_new, rho_original, 1, 'rho', 0.0001, 1, t, dt)
##  fl.Porta_plotty2(x, rho_newUpwind, rho_original, 1, 'rho', 0.0001, 2, t, dt)
#
  print(len(q))
##  plt.show()
#
##  fig, axs = plt.subplots(2)
#  if (MOVIE == 1):
#    plt.figure(1)
#    plt.clf()
#    plt.title(' Step number  = ' + str(t) )
#    plt.plot(x[:], rho_newUpwind, label = 'Density')
###    plt.plot(x[:], mom_new, label = 'Momentum')
##    plt.plot(x[:], vel, label='Velocity')
#    plt.pause(0.1)
#    plt.legend()
#    plt.figure(2)
#    plt.clf()
#    plt.title('Pressure--Step number == '+ str(t))
#    plt.plot(x[:], Pressure, label = 'Pressure')
##    plt.plot(x[:], rho_newUpwind, label = 'Upwind')
##    plt.plot(x[:], rho_newFromms, label = 'Fromms')
##    plt.plot(x[:], rho_original[1:-1])
##    plt.axhline(y=1, color='r', linestyle='--')
##    plt.axis((0, 2*np.pi, -1, 1))
#    plt.pause(0.5)
#    plt.legend()

#  if(t*dt > 71.):
#    plt.show()
#
#    plt.figure(2)
#    plt.clf()
#    plt.title('Momentum Time $t$ = ' + str(t*dt) + ' ns')
#    plt.plot(x[1:-1], mom_new) #, label = str((n-1) * 0.1) + ' ns')
##    plt.axhline(y=1, color='r', linestyle='--')
#    plt.axis((0, 2*np.pi, -1, 1))
#    plt.pause(0.0001)
#    plt.legend()
#  
#    plt.figure(3)
#    plt.clf()
#    plt.title('Energy Density Time $t$ = ' + str(t*dt) + ' ns')
#    plt.plot(x[1:-1], engyDens_new) #, label = str((n-1) * 0.1) + ' ns')
##    plt.axhline(y=1, color='r', linestyle='--')
##    plt.axis((0, 2*np.pi, -1, 1))
#    plt.pause(0.0001)
#    plt.legend()
  
#    plt.figure(4)
#    plt.clf()
#    plt.title('Temperature Time $t$ = ' + str(t*dt) + ' ns')
#    plt.plot(x[1:-1], temp) #, label = str((n-1) * 0.1) + ' ns')
##    plt.axhline(y=1, color='r', linestyle='--')
##    plt.axis((0, 2*np.pi))
#    plt.pause(0.0001)
#    plt.legend()
#
#    plt.figure(5)
#    plt.clf()
#    plt.title('Pressure Time $t$ = ' + str(t*dt) + ' ns')
#    plt.plot(x[1:], Pressure) #, label = str((n-1) * 0.1) + ' ns')
##    plt.axhline(y=1, color='r', linestyle='--')
##    plt.axis((0, 2*np.pi))
#    plt.pause(0.0001)
#    plt.legend()
#
#    plt.figure(6)
#    plt.clf()
#    plt.title('Velocity Time $t$ = ' + str(t*dt) + ' ns')
#    plt.plot(x[1:-1], vel) #, label = str((n-1) * 0.1) + ' ns')
##    plt.axhline(y=1, color='r', linestyle='--')
##    plt.axis((0, 2*np.pi))
#    plt.pause(0.0001)
#    plt.legend()

