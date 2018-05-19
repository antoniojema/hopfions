import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pylab
import h5py as h5

#Open data file
fin = h5.File('simulation_results.h5','r')

N = fin.attrs['N'][0]
iterations = fin.attrs['iterations'][0]
L = 10
x = np.linspace(0,L,N+1)
x_ = range(N+1)

def ind(i,j,k):
	return (N+1)*(N+1)*i+(N+1)*j+k

### Animation ###
fig = plt.figure()

def iteration(i):
	global x_
	global x
	plt.cla()
	plt.ylim([0,3])
	plt.xlim([0,10])
	E = np.sqrt((fin[str(i)]['Ex'][:])**2+(fin[str(i)]['Ey'][:])**2+(fin[str(i)]['Ez'][:])**2)
	plt.plot(x,E[ind(N/2,N/2,x_)])
	print(i)
	return

animation = ani.FuncAnimation(fig, iteration, iterations, interval=25)
plt.show()
