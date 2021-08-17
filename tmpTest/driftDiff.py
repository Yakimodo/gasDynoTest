import numpy as np
import matplotlib.pyplot as plt
import function_list_Test as fl



N = 1001
L = 1.
CFL = 0.6
dx = L / (2.*(N-1)) #spatial resolution
dt = CFL*dx #time step
T_begin = 1 #step number [important in determining xi = x/t ----> (in this case) --> xc/(t*dt)]
T_end = 500 #number of steps
gamma = 1.4 #gravity
xc = np.linspace(-L/2, L/2, N)

eqNum = 2

n_empty = np.empty((eqNum, N+2))

n = initCond('Gauss', xc, N, n_empty)
#n[0,:] = ne
#n[1,:] = ni
#n[2,:] = N (neutrals O2 and N2)
#
























