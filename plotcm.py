import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pylab
import h5py as h5

#Open data file
fin = h5.File('simulation_results.h5','r')

N = fin.attrs['N'][0]
iterations = fin.attrs['iterations'][0]
I = np.array([[i for j in range(N+1)] for i in range(N+1)])
J = np.array([[j for j in range(N+1)] for i in range(N+1)])

def ind(i,j,k):
	return (N+1)*(N+1)*i+(N+1)*j+k

### Animation ###
fig = plt.figure()

def init():
	global cbar
	
	#Ex = fin['0']['Ex'][:]
	#Ey = fin['0']['Ey'][:]
	Ez = fin['0']['Ez'][:]
	#Ex = Ex[ind(I,J,N/2)]	
	#Ey = Ey[ind(I,J,N/2)]
	#Ez = 0.5*(Ez[ind(I,J,N/2)]+Ez[ind(I,J,N/2+1)])
	#Ex = np.array([[0.5*(Ex[i][j] + Ex[i+1][j]) for j in range(N)] for i in range(N)])
	#Ey = np.array([[0.5*(Ey[i][j] + Ey[i][j+1]) for j in range(N)] for i in range(N)])
	#Ez = np.array([[Ez[i][j] for j in range(N)] for i in range(N)])
	#P = Ex*Ex+Ey*Ey+Ez*Ez
	
	#pylab.pcolor(P)
	Ez = Ez[ind(I,J,N/2)]
	pylab.pcolor(Ez)
	#plt.clim(0,1)
	cbar = pylab.colorbar()

def iteration(n):
	global cbar
	global I
	global J
	plt.cla()
	
	#Ex = fin[str(n)]['Ex'][:]
	#Ey = fin[str(n)]['Ey'][:]
	Ez = fin[str(n)]['Ez'][:]
	#Ex = Ex[ind(I,J,N/2)]	
	#Ey = Ey[ind(I,J,N/2)]
	#Ez = 0.5*(Ez[ind(I,J,N/2)]+Ez[ind(I,J,N/2+1)])
	#Ex = np.array([[0.5*(Ex[i][j] + Ex[i+1][j]) for j in range(N)] for i in range(N)])
	#Ey = np.array([[0.5*(Ey[i][j] + Ey[i][j+1]) for j in range(N)] for i in range(N)])
	#Ez = np.array([[Ez[i][j] for j in range(N)] for i in range(N)])
	#P = Ex*Ex+Ey*Ey+Ez*Ez
	
	print n
	#pylab.pcolor(P)
	Ez = Ez[ind(I,J,N/2)]
	pylab.pcolor(abs(Ez))
	#plt.clim(0,1)
	cbar.remove()
	cbar = pylab.colorbar()

animation = ani.FuncAnimation(fig, iteration, np.arange(5,iterations,5), interval=25, init_func=init)
plt.show()

