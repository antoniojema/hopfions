import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pylab
import h5py as h5
import time
from F23 import *

#Open data file
fin = h5.File('simulation_results.h5','r')

N = fin.attrs['N'][0]
iterations = fin.attrs['iterations'][0]
Dt = fin.attrs['Dt'][0]
D = 2.*Dt
I_z = np.array([[i for j in range(N+1)] for i in range(N+1)])
J_z = np.array([[j for j in range(N+1)] for i in range(N+1)])
I_P = np.array([[i for j in range(N)] for i in range(N)])
J_P = np.array([[j for j in range(N)] for i in range(N)])
I = D * (I_z - 0.5*N)
J = D * (J_z - 0.5*N)

def ind(i,j,k):
	global N
	return (N+1)*(N+1)*i+(N+1)*j+k

### Animation ###
fig, ax = plt.subplots()

def init():
	global cbar, I_P, J_P, I_z, I_z, I, J, ax
	print 0
	
	Ex = fin['0']['Ex'][:]
	Ey = fin['0']['Ey'][:]
	Ez = fin['0']['Ez'][:]
	
	Ex = 0.5 * ( Ex[ind(I_P,N/2,J_P)] + Ex[ind(I_P+1, N/2 , J_P )] )
	Ey = 0.5 * ( Ey[ind(I_P,N/2,J_P)] + Ey[ind( I_P ,N/2+1, J_P )] )
	Ez = 0.5 * ( Ez[ind(I_P,N/2,J_P)] + Ez[ind( I_P , N/2 ,J_P+1)] )
	
	P = Ex*Ex + Ey*Ey + Ez*Ez
	
	pylab.pcolor(P)
	plt.clim(0,20)
	
	'''
	Ez = fin['0']['Ez'][:]
	Ez = Ez[ind(I_z,J_z,N/2)]
	
	pylab.pcolor(abs(Ez))
	plt.clim(0,5)
	'''
	'''
	Ex = ( Fx(    I    , 0.5*D , J+0.5*D , -1.5-0.5*Dt ) ).real
	Ey = ( Fy( I+0.5*D ,   0   , J+0.5*D , -1.5-0.5*Dt ) ).real
	Ez = ( Fz( I+0.5*D , 0.5*D ,    J    , -1.5-0.5*Dt ) ).real
	
	P = Ex*Ex + Ey*Ey + Ez*Ez
	
	pylab.pcolor(P)
	plt.clim(0,20)
	'''
	'''
	Ex = ( Fx( I+0.5*D , 0.5*D , J+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
	
	pylab.pcolor(abs(Ex))
	plt.clim(0,5)
	'''
	cbar = pylab.colorbar()
	
	ax.text(1,190,'t = '+"%.3f"%(-1.5-0.5*Dt + 0*Dt)+' s',fontsize=12,color='white')

def iteration(n):
	global cbar, I_P, J_P, I_z, I_z, ax
	plt.cla()
	print n
	
	Ex = fin[str(n)]['Ex'][:]
	Ey = fin[str(n)]['Ey'][:]
	Ez = fin[str(n)]['Ez'][:]
	
	Ex = 0.5 * ( Ex[ind(I_P,N/2,J_P)] + Ex[ind(I_P+1, N/2 , J_P )] )
	Ey = 0.5 * ( Ey[ind(I_P,N/2,J_P)] + Ey[ind( I_P ,N/2+1, J_P )] )
	Ez = 0.5 * ( Ez[ind(I_P,N/2,J_P)] + Ez[ind( I_P , N/2 ,J_P+1)] )
	
	P = Ex*Ex + Ey*Ey + Ez*Ez
	
	pylab.pcolor(P)
	plt.clim(0,20)
	
	'''
	Ez = fin[str(n)]['Ez'][:]
	Ez = Ez[ind(I_z,J_z,N/2)]
	
	pylab.pcolor(abs(Ez))
	plt.clim(0,5)
	'''
	'''
	Ex = ( Fx(    I    , 0.5*D , J+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
	Ey = ( Fy( I+0.5*D ,   0   , J+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
	Ez = ( Fz( I+0.5*D , 0.5*D ,    J    , -1.5-0.5*Dt+n*Dt ) ).real
	
	P = Ex*Ex + Ey*Ey + Ez*Ez
	
	pylab.pcolor(P)
	plt.clim(0,20)
	'''
	'''
	Ex = ( Fx( I+0.5*D , 0.5*D , J+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
	
	pylab.pcolor(abs(Ex))
	plt.clim(0,5)
	'''
	cbar.remove()
	cbar = pylab.colorbar()
	
	ax.text(1,190,'t = '+"%.3f"%(-1.5-0.5*Dt + n*Dt)+' s',fontsize=12,color='white')

animation = ani.FuncAnimation(fig, iteration, np.arange(1,iterations+1,1), interval=25, init_func=init)
#plt.show()

t0 = time.time()
Writer = ani.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
animation.save('videos/P_XZ_simulac.mp4',writer=writer)
print time.time()-t0,' s'

