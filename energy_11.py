import numpy as np
import h5py as h5
from mayavi import mlab
from F11 import *
import sys

SEEDS = 8
n=int(sys.argv[2])
ARG1 = sys.argv[1]

def ind(i,j,k):
	global N
	return (N+1)*(N+1)*i+(N+1)*j+k


if ARG1 == 'sim':
	###### INPUT DATA ######
	
	fin = h5.File('results_11.h5','r')
	N = fin.attrs['N'][0]
	Dt = fin.attrs['Dt'][0]
	
	I = np.array([[[i for k in range(N)] for j in range(N)] for i in range(N)])+1
	J = np.array([[[j for k in range(N)] for j in range(N)] for i in range(N)])+1
	K = np.array([[[k for k in range(N)] for j in range(N)] for i in range(N)])+1
	
	mlab.figure(size=(1366,1366),fgcolor=(0,0,0),bgcolor=(1,1,1))
	
	Ex = fin[str(n)]['Ex'][:]
	Ey = fin[str(n)]['Ey'][:]
	Ez = fin[str(n)]['Ez'][:]
	
	Ex = 0.25  *( Ex[ind(I,J,K)] + Ex[ind( I ,J-1, K )] + Ex[ind( I , J ,K-1)] + Ex[ind( I ,J-1,K-1)] )
	Ey = 0.25  *( Ey[ind(I,J,K)] + Ey[ind(I-1, J , K )] + Ey[ind( I , J ,K-1)] + Ey[ind(I-1, J ,K-1)] )
	Ez = 0.25  *( Ez[ind(I,J,K)] + Ez[ind(I-1, J , K )] + Ez[ind( I ,J-1, K )] + Ez[ind(I-1,J-1, K )] )
	
	Hx = fin[str(n)]['Hx'][:]
	Hy = fin[str(n)]['Hy'][:]
	Hz = fin[str(n)]['Hz'][:]
	
	Hx1 = fin[str(n-1)]['Hx'][:]
	Hy1 = fin[str(n-1)]['Hy'][:]
	Hz1 = fin[str(n-1)]['Hz'][:]
	
	Hx = 0.25  *( Hx[ind(I,J,K)] + Hx[ind(I-1, J , K )] + Hx1[ind(I,J,K)] + Hx1[ind(I-1, J , K )] )
	Hy = 0.25  *( Hy[ind(I,J,K)] + Hy[ind( I ,J-1, K )] + Hy1[ind(I,J,K)] + Hy1[ind( I ,J-1, K )] )
	Hz = 0.25  *( Hz[ind(I,J,K)] + Hz[ind( I , J ,K-1)] + Hz1[ind(I,J,K)] + Hz1[ind( I , J ,K-1)] )
	del Hx1, Hy1, Hz1
	
	Px = Ey*Hz-Ez*Hy
	Py = Ez*Hx-Ex*Hz
	Pz = Ex*Hy-Ey*Hx
	del Ex, Ey, Ez, Hx, Hy, Hz
	
	P = np.sqrt(Px*Px+Py*Py+Pz*Pz)
	module = mlab.pipeline.scalar_field(P); del Px, Py, Pz, P;
	mlab.text(0.65,0.9,'t = '+"%.4f"%(-1.5-0.5*Dt + n*Dt)+' s',width=0.3)
	color = (0,1,0)
	
	mlab.pipeline.volume(module, vmin=0, vmax=0.8)
	
	mlab.outline(extent=[0,200,0,200,0,200],opacity=0.2,color=(0,0,0))
	mlab.orientation_axes()
	
	
	'''
	mlab.view(distance=700)
	if n == 1:
		mlab.savefig('../hopfions_memoria/media/11_energy_t-15_'+ARG1+'.png',size=(1366,1366))
	elif n == 60:
		mlab.savefig('../hopfions_memoria/media/11_energy_t0_'+ARG1+'.png',size=(1366,1366))
	elif n == 120:
		mlab.savefig('../hopfions_memoria/media/11_energy_t15_'+ARG1+'.png',size=(1366,1366))
	'''
	
	#mlab.savefig('../hopfions_memoria/video/11_energy_'+str(n)+'_'+ARG1+'.png',size=(1366,1366))
	
	mlab.show()

###################################################################### ^ SIMULADO ^


#                        ~Territorio de nadie~						 #


###################################################################### v TEORICO  v

elif ARG1 == 'teor':
	###### INPUT DATA ######
	
	N=200
	D = 10./N
	Dt = 0.5*D
	
	I = D * (np.array([[[i for k in range(N+1)] for j in range(N+1)] for i in range(N+1)]) - 0.5*N)
	J = D * (np.array([[[j for k in range(N+1)] for j in range(N+1)] for i in range(N+1)]) - 0.5*N)
	K = D * (np.array([[[k for k in range(N+1)] for j in range(N+1)] for i in range(N+1)]) - 0.5*N)
	
	mlab.figure(size=(1366,1366),fgcolor=(0,0,0),bgcolor=(1,1,1))
	
	Ex = ( Fx( I , J , K , -1.5-0.5*Dt+n*Dt ) ).real
	Ey = ( Fy( I , J , K , -1.5-0.5*Dt+n*Dt ) ).real
	Ez = ( Fz( I , J , K , -1.5-0.5*Dt+n*Dt ) ).real
	
	Hx = ( Fx( I , J , K , -1.5-0.5*Dt+n*Dt ) ).imag
	Hy = ( Fy( I , J , K , -1.5-0.5*Dt+n*Dt ) ).imag
	Hz = ( Fz( I , J , K , -1.5-0.5*Dt+n*Dt ) ).imag
	
	Px = Ey*Hz-Ez*Hy
	Py = Ez*Hx-Ex*Hz
	Pz = Ex*Hy-Ey*Hx
	del Ex, Ey, Ez, Hx, Hy, Hz
	
	P = np.sqrt(Px*Px+Py*Py+Pz*Pz)
	module = mlab.pipeline.scalar_field(P); del Px, Py, Pz, P;
	mlab.text(0.65,0.9,'t = '+"%.4f"%(-1.5-0.5*Dt + n*Dt)+' s',width=0.3)
	color = (0,1,0)
	
	mlab.pipeline.volume(module, vmin=0, vmax=0.8)
	
	mlab.outline(extent=[0,200,0,200,0,200],opacity=0.2,color=(0,0,0))
	mlab.orientation_axes()
	
	
	
	'''
	mlab.view(distance=700)
	if n == 1:
		mlab.savefig('../hopfions_memoria/media/11_energy_t-15_'+ARG1+'.png',size=(1366,1366))
	elif n == 60:
		mlab.savefig('../hopfions_memoria/media/11_energy_t0_'+ARG1+'.png',size=(1366,1366))
	elif n == 120:
		mlab.savefig('../hopfions_memoria/media/11_energy_t15_'+ARG1+'.png',size=(1366,1366))
	'''
	
	#mlab.savefig('../hopfions_memoria/video/11_energy_'+str(n)+'_'+ARG1+'.png',size=(1366,1366))
	
	mlab.show()

else:
	print 'ARGUMENT ERROR'