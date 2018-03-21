import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pylab
import h5py as h5

#Open data file
fin = h5.File('3D.h5','r')

N=100
iterations = 200
'''
N = fin.attrs['N']
iterations = fin.attrs['iterations']
'''
L = 10
x = np.linspace(0,L,N+1)

### Animation ###
fig = plt.figure()
plt.ylim([0,3])
plt.xlim([0,10])

def iteration(i):
	plt.cla();
	plt.autoscale(False)
	plt.plot(x,np.sqrt((fin[str(i)]['Ex'][N/2][N/2][:])**2+(fin[str(i)]['Ey'][N/2][N/2][:])**2+(fin[str(i)]['Ez'][N/2][N/2][:])**2))
	return

animation = ani.FuncAnimation(fig, iteration, iterations, interval=25)
plt.show();
