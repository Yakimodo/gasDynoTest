import matplotlib.pyplot as plt
import numpy as np
import csv
import function_list as fl

#What is needed? 
# Solution to mass continuity equation --- 1D Cartesian : dp/dt = - div ( rho * vel ) 
# ICs = Density and Velocity
MOVIE = 1 # 1=yes 0=no

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

##Constants----
N = 100 #number of grid cells
ng = 2 #number of ghost cells
CFL = 0.00001 #0.0001
dx = 1./N
dt = CFL*dx
beta = 200.
T_end = 200
x = np.linspace(0, 2*np.pi, N+2)
xmax = max(x)
xmin = min(x)
kB = 1.38e-23
Ndens = 1e21


##----Evolve in time----
for t in range(0,T_end):
  print(t)



##---Initialize-----
  if (t == 0):
    rho, vel, engyDens, Pressure, rho_original = fl.InitialConditions(x,N)


##--Interface states----
  Qr_cont = []
  Ql_cont = []
  Qr_mom = []
  Ql_mom = []
  Qr_engyDens = []
  Ql_engyDens = []
  Qt = []

  rho = fl.fill_ghosts(rho, ng, N)
  vel = fl.fill_ghosts(vel, ng, N)
  engyDens = fl.fill_ghosts(engyDens, ng, N)
  rho1 = rho
  Qt = fl.PowerDeposition(N, dt)

###---initializing temperature-------###
  if (t == 0):
    temp = fl.Temp(Ndens, engyDens, rho, vel, N) #initialize temperature
  temp = fl.fill_ghosts(temp, ng, N)

###----intializaing momentum------###
  for i in range(len(vel)):
    mom.append(rho[i] * vel[i])

  gradP = fl.Upwind(Pressure, dx, N)
  Qr_cont, Ql_cont = fl.states(rho, dx, dt, vel, 'SuperBee', N)
  Qr_mom, Ql_mom = fl.states(mom, dx, dt, vel, 'SuperBee', N)
  Qr_engyDens, Ql_engyDens = fl.states(engyDens, dx, dt, vel, 'SuperBee', N)
  
##--Solve Jump Discontinuity----
  Flux_cont = []
  Flux_mom = []
  Flux_engyDens = []
  Flux_cont = fl.Riemann(Ql_cont, Qr_cont, vel, N)
  Flux_mom = fl.Riemann(Ql_mom, Qr_mom, vel, N)
  Flux_engyDens = fl.Riemann(Ql_engyDens, Qr_engyDens, vel, N)
 
##--To update solution----
  rho_new = []
  rho_new1 = []
  mom_new = [] 
  engyDens_new = []
 
  for i in range(N-1):
     rho_new.append(rho[i] - (dt/dx) * (Flux_cont[i] - Flux_cont[i+1]))
     mom_new.append(mom[i] - (dt/dx) * (Flux_mom[i] - Flux_mom[i+1]) - dt * gradP[i])
     rho_new1.append(rho1[i] * (1 - 2*(dt/dx)*vel[i]) + (dt/dx)*(vel[i]*rho1[i-1] + rho1[i]*vel[i-1]))
     engyDens_new.append(engyDens[i] - (dt/dx)*(Flux_engyDens[i] - Flux_engyDens[i+1]) + dt*Qt)

##---accounting for ghosts---##
  rho_new.insert(0, rho_new[N-2])
  rho_new1.insert(0, rho_new[N-2])
  rho = rho_new
  rho1 = rho_new1
  engyDens_new.insert(0, engyDens_new[N-2])
  mom_new.insert(0, mom_new[N-2])
  mom = mom_new

##---updating temperature-----##
  temp = fl.Temp(N, engyDens, rho, vel)

##---updating pressure------##
  Pressure = []
  for i in range(N-1):
    Pressure.append(Ndens * kB * temp[i])


##---velocity matters----##
  vel = vel[1:-1]

###-----PLOT-----#####
#  fig, axs = plt.subplots(2)
  if (MOVIE == 1):
    plt.figure(1)
    plt.clf()
    plt.title('Time $t$ = ' + str(t*dt) + ' ns')
    plt.plot(x[1:-1], rho_new) #, label = str((n-1) * 0.1) + ' ns')
    plt.plot(x[1:-1], rho_new1) #, label = str((n-1) * 0.1) + ' ns')
    plt.plot(x, rho_original[:])
#    plt.axhline(y=1, color='r', linestyle='--')
    plt.axis((0, 2*np.pi, -1, 1))
    plt.pause(0.0001)
    plt.legend()

    plt.figure(2)
    plt.clf()
    plt.title('Time $t$ = ' + str(t*dt) + ' ns')
    plt.plot(x[1:-1], mom_new) #, label = str((n-1) * 0.1) + ' ns')
#    plt.axhline(y=1, color='r', linestyle='--')
    plt.axis((0, 2*np.pi, -1, 1))
    plt.pause(0.0001)
    plt.legend()
  
    plt.figure(3)
    plt.clf()
    plt.title('Time $t$ = ' + str(t*dt) + ' ns')
    plt.plot(x[1:-1], engyDens_new) #, label = str((n-1) * 0.1) + ' ns')
#    plt.axhline(y=1, color='r', linestyle='--')
    plt.axis((0, 2*np.pi, -1, 1))
    plt.pause(0.0001)
    plt.legend()
  



