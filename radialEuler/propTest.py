	
import numpy as np
import matplotlib.pyplot as plt
import function_list_Test as fl




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


EfOld = []
for j in range(N):
  EfOld.append(np.exp(-np.abs((rc[j]-rc[(N-1)/4]))**2 / 0.02))
#EfOld = fl.fill_ghosts(EfOld, N, 'extrap')
#EfOld = np.array(EfOld)

JOld = []
for j in range(N):
  JOld.append(np.exp(-np.abs((rc[j]-rc[np.int(3.*(N-1)//4)]))**2 / 0.02))


for t in range(T_end):
  EfOld = fl.fill_ghosts(EfOld, N, 'extrap')
  JOld = fl.fill_ghosts(JOld, N, 'extrap')

  c = 1.
  EfNew = fl.simplePropRight(EfOld, dt, dr, 1.)
  JNew = fl.simplePropLeft(JOld, dt, dr, -1.)

  plt.figure(2)
  plt.clf()
  plt.plot(rc, EfOld[:len(rc)], label = 'EF')
  plt.plot(rc, JOld[:len(rc)], label = 'Current Density')
#      plt.plot(rc, Temp[1:], label = 'Temperature')i
  plt.pause(0.1)
  plt.legend()

  EfOld = EfNew
  JOld = JNew



