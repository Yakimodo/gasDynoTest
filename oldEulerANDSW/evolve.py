import numpy as np
import matplotlib.pyplot as plt

DIR = '/home/yakimoto/Code/GasDyno_Test/'
#sim_type = 'Sound'
#sim_details = 'VImp'

sim_type = 'ShallowWater' #'Sound'
sim_details = 'SmallDisturbance' #'VImp'

arr1 = []
arr2 = []
N = 501
xc = np.linspace(0,10,N)

arr1 = np.load(DIR + 'Qevolv1' + sim_type + sim_details + '.npy')
arr2 = np.load(DIR + 'Qevolv2' + sim_type + sim_details + '.npy')

for n in range(len(arr1)):
  arr1 = np.array(arr1)
  arr2 = np.array(arr2)
  
  plt.figure(1)
  plt.clf()
  plt.title(' Step number  = ' + str(n) )
  plt.plot(xc, arr1[n,1:len(xc)+1])
  plt.plot(xc, arr2[n,1:len(xc)+1])
  if (sim_type == 'Sound'):
    plt.axvline(x=0.5*max(xc))
    plt.axis([0, 10, -1., 2]) 
  if (sim_type == 'ShallowWater'):
    plt.axis([0, 10, -1., 2]) 
  plt.pause(0.1)

