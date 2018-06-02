import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pylab
import h5py as h5
import time
from F23 import *
import sys

field = sys.argv[1]
plane = sys.argv[2]

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
I = D * (I_P - 0.5*N)
J = D * (J_P - 0.5*N)
I_ = D * (I_z - 0.5*N)
J_ = D * (J_z - 0.5*N)

def ind(i,j,k):
	global N
	return (N+1)*(N+1)*i+(N+1)*j+k

### Animation ###
fig, ax = plt.subplots()

def init():
	global cbar, I_P, J_P, I_z, I_z, I, J, ax, field, plane
	#print 0
	
	if field=='P_sim':
		Ex = fin['0']['Ex'][:]
		Ey = fin['0']['Ey'][:]
		Ez = fin['0']['Ez'][:]
		if plane=='XY':
			Ex = 0.5 * ( Ex[ind(I_P,J_P,N/2)] + Ex[ind(I_P+1, J_P , N/2 )] )
			Ey = 0.5 * ( Ey[ind(I_P,J_P,N/2)] + Ey[ind( I_P ,J_P+1, N/2 )] )
			Ez = 0.5 * ( Ez[ind(I_P,J_P,N/2)] + Ez[ind( I_P , J_P ,N/2+1)] )
		elif plane=='XZ':
			Ex = 0.5 * ( Ex[ind(I_P,N/2,J_P)] + Ex[ind(I_P+1, N/2 , J_P )] )
			Ey = 0.5 * ( Ey[ind(I_P,N/2,J_P)] + Ey[ind( I_P ,N/2+1, J_P )] )
			Ez = 0.5 * ( Ez[ind(I_P,N/2,J_P)] + Ez[ind( I_P , N/2 ,J_P+1)] )
		elif plane=='YZ':
			Ex = 0.5 * ( Ex[ind(N/2,I_P,J_P)] + Ex[ind(N/2+1, I_P , J_P )] )
			Ey = 0.5 * ( Ey[ind(N/2,I_P,J_P)] + Ey[ind( N/2 ,I_P+1, J_P )] )
			Ez = 0.5 * ( Ez[ind(N/2,I_P,J_P)] + Ez[ind( N/2 , I_P ,J_P+1)] )
	
		P = Ex*Ex + Ey*Ey + Ez*Ez
		pylab.pcolor(P)
		plt.clim(0,20)
	
	elif field=='P_teor':
		if plane=='XY':
			Ex = ( Fx( I+0.5*D , J+0.5*D , 0+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
			Ey = ( Fy( I+0.5*D , J+0.5*D , 0+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
			Ez = ( Fz( I+0.5*D , J+0.5*D , 0+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		elif plane=='XZ':
			Ex = ( Fx( I+0.5*D , 0+0.5*D , J+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
			Ey = ( Fy( I+0.5*D , 0+0.5*D , J+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
			Ez = ( Fz( I+0.5*D , 0+0.5*D , J+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		elif plane=='YZ':
			Ex = ( Fx( 0+0.5*D , I+0.5*D , J+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
			Ey = ( Fy( 0+0.5*D , I+0.5*D , J+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
			Ez = ( Fz( 0+0.5*D , I+0.5*D , J+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		
		P = Ex*Ex + Ey*Ey + Ez*Ez
		pylab.pcolor(P)
		plt.clim(0,20)
	
	elif field=='Ex_sim':
		E = fin['0']['Ex'][:]
		if plane=='XY':
			E = E[ind(I_z,J_z,N/2)]
		elif plane=='XZ':
			E = E[ind(I_z,N/2,J_z)]
		elif plane=='YZ':
			E = E[ind(N/2,I_z,J_z)]
		
		pylab.pcolor(abs(E))
		plt.clim(0,5)
	
	elif field =='Ex_teor':
		if plane=='XY':
			E = ( Fx( I_ , J_+0.5*D ,  0+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		elif plane=='XZ':
			E = ( Fx( I_ ,  0+0.5*D , J_+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		elif plane=='YZ':
			E = ( Fx(  0 , I_+0.5*D , J_+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		
		pylab.pcolor(abs(E))
		plt.clim(0,5)
	
	elif field=='Ey_sim':
		E = fin['0']['Ey'][:]
		if plane=='XY':
			E = E[ind(I_z,J_z,N/2)]
		elif plane=='XZ':
			E = E[ind(I_z,N/2,J_z)]
		elif plane=='YZ':
			E = E[ind(N/2,I_z,J_z)]
		
		pylab.pcolor(abs(E))
		plt.clim(0,5)
	
	elif field =='Ey_teor':
		if plane=='XY':
			E = ( Fy( I_+0.5*D , J_ ,  0+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		elif plane=='XZ':
			E = ( Fy( I_+0.5*D ,  0 , J_+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		elif plane=='YZ':
			E = ( Fy(  0+0.5*D , I_ , J_+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		
		pylab.pcolor(abs(E))
		plt.clim(0,5)
	
	elif field=='Ez_sim':
		E = fin['0']['Ez'][:]
		if plane=='XY':
			E = E[ind(I_z,J_z,N/2)]
		elif plane=='XZ':
			E = E[ind(I_z,N/2,J_z)]
		elif plane=='YZ':
			E = E[ind(N/2,I_z,J_z)]
		
		pylab.pcolor(abs(E))
		plt.clim(0,5)
	
	elif field =='Ez_teor':
		if plane=='XY':
			E = ( Fz( I_+0.5*D , J_+0.5*D ,  0 , -1.5-0.5*Dt+0*Dt ) ).real
		elif plane=='XZ':
			E = ( Fz( I_+0.5*D ,  0+0.5*D , J_ , -1.5-0.5*Dt+0*Dt ) ).real
		elif plane=='YZ':
			E = ( Fz(  0+0.5*D , I_+0.5*D , J_ , -1.5-0.5*Dt+0*Dt ) ).real
		
		pylab.pcolor(abs(E))
		plt.clim(0,5)
	
	cbar = pylab.colorbar()
	
	ax.text(1,190,'t = '+"%.3f"%(-1.5-0.5*Dt + 0*Dt)+' s',fontsize=12,color='white')

def iteration(n):
	global cbar, I_P, J_P, I_z, I_z, ax, field, plane
	plt.cla()
	#print n
	
	if field=='P_sim':
		Ex = fin[str(n)]['Ex'][:]
		Ey = fin[str(n)]['Ey'][:]
		Ez = fin[str(n)]['Ez'][:]
		if plane=='XY':
			Ex = 0.5 * ( Ex[ind(I_P,J_P,N/2)] + Ex[ind(I_P+1, J_P , N/2 )] )
			Ey = 0.5 * ( Ey[ind(I_P,J_P,N/2)] + Ey[ind( I_P ,J_P+1, N/2 )] )
			Ez = 0.5 * ( Ez[ind(I_P,J_P,N/2)] + Ez[ind( I_P , J_P ,N/2+1)] )
		elif plane=='XZ':
			Ex = 0.5 * ( Ex[ind(I_P,N/2,J_P)] + Ex[ind(I_P+1, N/2 , J_P )] )
			Ey = 0.5 * ( Ey[ind(I_P,N/2,J_P)] + Ey[ind( I_P ,N/2+1, J_P )] )
			Ez = 0.5 * ( Ez[ind(I_P,N/2,J_P)] + Ez[ind( I_P , N/2 ,J_P+1)] )
		elif plane=='YZ':
			Ex = 0.5 * ( Ex[ind(N/2,I_P,J_P)] + Ex[ind(N/2+1, I_P , J_P )] )
			Ey = 0.5 * ( Ey[ind(N/2,I_P,J_P)] + Ey[ind( N/2 ,I_P+1, J_P )] )
			Ez = 0.5 * ( Ez[ind(N/2,I_P,J_P)] + Ez[ind( N/2 , I_P ,J_P+1)] )
	
		P = Ex*Ex + Ey*Ey + Ez*Ez
		pylab.pcolor(P)
		plt.clim(0,20)
	
	elif field=='P_teor':
		if plane=='XY':
			Ex = ( Fx( I+0.5*D , J+0.5*D , 0+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
			Ey = ( Fy( I+0.5*D , J+0.5*D , 0+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
			Ez = ( Fz( I+0.5*D , J+0.5*D , 0+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		elif plane=='XZ':
			Ex = ( Fx( I+0.5*D , 0+0.5*D , J+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
			Ey = ( Fy( I+0.5*D , 0+0.5*D , J+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
			Ez = ( Fz( I+0.5*D , 0+0.5*D , J+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		elif plane=='YZ':
			Ex = ( Fx( 0+0.5*D , I+0.5*D , J+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
			Ey = ( Fy( 0+0.5*D , I+0.5*D , J+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
			Ez = ( Fz( 0+0.5*D , I+0.5*D , J+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		
		P = Ex*Ex + Ey*Ey + Ez*Ez
		pylab.pcolor(P)
		plt.clim(0,20)
	
	elif field=='Ex_sim':
		E = fin[str(n)]['Ex'][:]
		if plane=='XY':
			E = E[ind(I_z,J_z,N/2)]
		elif plane=='XZ':
			E = E[ind(I_z,N/2,J_z)]
		elif plane=='YZ':
			E = E[ind(N/2,I_z,J_z)]
		
		pylab.pcolor(abs(E))
		plt.clim(0,5)
	
	elif field =='Ex_teor':
		if plane=='XY':
			E = ( Fx( I_ , J_+0.5*D ,  0+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		elif plane=='XZ':
			E = ( Fx( I_ ,  0+0.5*D , J_+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		elif plane=='YZ':
			E = ( Fx(  0 , I_+0.5*D , J_+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		
		pylab.pcolor(abs(E))
		plt.clim(0,5)
	
	elif field=='Ey_sim':
		E = fin[str(n)]['Ey'][:]
		if plane=='XY':
			E = E[ind(I_z,J_z,N/2)]
		elif plane=='XZ':
			E = E[ind(I_z,N/2,J_z)]
		elif plane=='YZ':
			E = E[ind(N/2,I_z,J_z)]
		
		pylab.pcolor(abs(E))
		plt.clim(0,5)
	
	elif field =='Ey_teor':
		if plane=='XY':
			E = ( Fy( I_+0.5*D , J_ ,  0+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		elif plane=='XZ':
			E = ( Fy( I_+0.5*D ,  0 , J_+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		elif plane=='YZ':
			E = ( Fy(  0+0.5*D , I_ , J_+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		
		pylab.pcolor(abs(E))
		plt.clim(0,5)
	
	elif field=='Ez_sim':
		E = fin[str(n)]['Ez'][:]
		if plane=='XY':
			E = E[ind(I_z,J_z,N/2)]
		elif plane=='XZ':
			E = E[ind(I_z,N/2,J_z)]
		elif plane=='YZ':
			E = E[ind(N/2,I_z,J_z)]
		
		pylab.pcolor(abs(E))
		plt.clim(0,5)
	
	elif field =='Ez_teor':
		if plane=='XY':
			E = ( Fz( I_+0.5*D , J_+0.5*D ,  0 , -1.5-0.5*Dt+n*Dt ) ).real
		elif plane=='XZ':
			E = ( Fz( I_+0.5*D ,  0+0.5*D , J_ , -1.5-0.5*Dt+n*Dt ) ).real
		elif plane=='YZ':
			E = ( Fz(  0+0.5*D , I_+0.5*D , J_ , -1.5-0.5*Dt+n*Dt ) ).real
		
		pylab.pcolor(abs(E))
		plt.clim(0,5)
	
	cbar.remove()
	cbar = pylab.colorbar()
	
	ax.text(1,190,'t = '+"%.3f"%(-1.5-0.5*Dt + n*Dt)+' s',fontsize=12,color='white')

animation = ani.FuncAnimation(fig, iteration, np.arange(1,iterations+1,1), interval=25, init_func=init)
#plt.show()

t0 = time.time()
Writer = ani.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
animation.save('videos/'+field+'_'+plane+'.mp4',writer=writer)
print time.time()-t0,' s'

