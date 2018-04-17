import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pylab
import h5py as h5

#Open data file
fin = h5.File('3D.h5','r')

N = fin.attrs['N']
iterations = fin.attrs['iterations']

def ind(i,j,k):
	return (N+1)*(N+1)*i+(N+1)*j+k

### Animation ###
fig = plt.figure()
plt.ylim([0,N])
plt.xlim([0,N])

y,x = np.meshgrid(range(N+1),range(N+1))

def init():
	global x
	global y
	global cbar
	E = np.sqrt((fin["0"]['Ex'][:])**2)
	pylab.pcolor(E[ind(x,y,N/2)])
	plt.clim(0,1)
	cbar = pylab.colorbar()

def iteration(i):
	global x
	global y
	global cbar
	plt.cla()
	
	plt.ylim([0,N])
	plt.xlim([0,N])
	E = np.absolute(fin[str(i)]['Ex'][:])
	pylab.pcolor(E[ind(x,y,50)])
	plt.clim(0,0.1)
	cbar.remove()
	cbar = pylab.colorbar()

animation = ani.FuncAnimation(fig, iteration, np.arange(1,iterations/5+1,1)*5, interval=25, init_func=init)
plt.show()

