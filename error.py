import numpy as np
import h5py as h5
from mayavi import mlab
from F11 import *
import sys

n=170
ARG = sys.argv[1]

def ind(i,j,k):
	global N
	return (N+1)*(N+1)*i+(N+1)*j+k

fin = h5.File('results_11.h5','r')
N = fin.attrs['N'][0]
D = fin.attrs['Dx'][0]
Dt = fin.attrs['Dt'][0]

I = np.array([[[i for k in range(N)] for j in range(N)] for i in range(N)])+1
J = np.array([[[j for k in range(N)] for j in range(N)] for i in range(N)])+1
K = np.array([[[k for k in range(N)] for j in range(N)] for i in range(N)])+1

if ARG == 'E':
	Ex = fin[str(n)]['Ex'][:]
	Ey = fin[str(n)]['Ey'][:]
	Ez = fin[str(n)]['Ez'][:]

	Ex = 0.5  *( Ex[ind(I-1,J-1,K-1)] + Ex[ind( I ,J-1,K-1)] )
	Ey = 0.5  *( Ey[ind(I-1,J-1,K-1)] + Ey[ind(I-1, J ,K-1)] )
	Ez = 0.5  *( Ez[ind(I-1,J-1,K-1)] + Ez[ind(I-1,J-1, K )] )

	Ex_ = ( Fx( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
	Ey_ = ( Fy( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
	Ez_ = ( Fz( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
	
	max_val = (np.sqrt(Ex_*Ex_+Ey_*Ey_+Ez_*Ez_)).max()
	E = np.sqrt( (Ex - Ex_)**2 + (Ey - Ey_)**2 + (Ez - Ez_)**2 )
	del Ex, Ey, Ez, Ex_, Ey_, Ez_
	
	mlab.figure(size=(500,500),fgcolor=(0,0,0),bgcolor=(1,1,1))
	mlab.pipeline.volume(mlab.pipeline.scalar_field(E), vmin=0.2, vmax=1)
	mlab.outline(extent=[0,200,0,200,0,200],opacity=0.2,color=(0,0,0))
	mlab.text(0.65,0.9,'t = '+'%.4f'%(-1.5-0.5*Dt + n*Dt)+' s',width=0.3)
	mlab.text(0.5,0.8,'Max. teor. = '+'%.2f'%max_val,width=0.45)


elif ARG == 'H':
	Hx = fin[str(n)]['Hx'][:]
	Hy = fin[str(n)]['Hy'][:]
	Hz = fin[str(n)]['Hz'][:]
	
	Hx = 0.5  *( Hx[ind(I,J,K)] + Hx[ind(I-1, J , K )] )
	Hy = 0.5  *( Hy[ind(I,J,K)] + Hy[ind( I ,J-1, K )] )
	Hz = 0.5  *( Hz[ind(I,J,K)] + Hz[ind( I , J ,K-1)] )
	
	Hx = ( Fx( I , J , K , -1.5+n*Dt ) ).imag
	Hy = ( Fy( I , J , K , -1.5+n*Dt ) ).imag
	Hz = ( Fz( I , J , K , -1.5+n*Dt ) ).imag
	
	max_val = (np.sqrt(Ex_*Ex_+Ey_*Ey_+Ez_*Ez_)).max()
	E = np.sqrt( (Hx - Hx_)**2 + (Hy - Hy_)**2 + (Hz - Hz_)**2 )
	del Hx, Hy, Hz, Hx_, Hy_, Hz_
	
	mlab.figure(size=(500,500),fgcolor=(0,0,0),bgcolor=(1,1,1))
	mlab.pipeline.volume(mlab.pipeline.scalar_field(E), vmin=0, vmax=1)
	mlab.outline(extent=[0,200,0,200,0,200],opacity=0.2,color=(0,0,0))
	mlab.text(0.65,0.9,'t = '+'%.4f'%(-1.5-0.5*Dt + n*Dt)+' s',width=0.3)
	mlab.text(0.5,0.8,'Max. teor. = '+'%.2f'%max_val,width=0.45)


elif ARG == 'P':
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
	
	Ex = ( Fx( I , J , K , -1.5-0.5*Dt+n*Dt ) ).real
	Ey = ( Fy( I , J , K , -1.5-0.5*Dt+n*Dt ) ).real
	Ez = ( Fz( I , J , K , -1.5-0.5*Dt+n*Dt ) ).real
	
	Hx = ( Fx( I , J , K , -1.5-0.5*Dt+n*Dt ) ).imag
	Hy = ( Fy( I , J , K , -1.5-0.5*Dt+n*Dt ) ).imag
	Hz = ( Fz( I , J , K , -1.5-0.5*Dt+n*Dt ) ).imag
	
	Px_ = Ey*Hz-Ez*Hy
	Py_ = Ez*Hx-Ex*Hz
	Pz_ = Ex*Hy-Ey*Hx
	del Ex, Ey, Ez
	
	max_val = (np.sqrt(Px_*Px_+Py_*Py_+Pz_*Pz_)).max()
	E = np.sqrt( (Px - Px_)**2 + (Py - Py_)**2 + (Pz - Pz_)**2 )
	del Px, Py, Pz, Px_, Py_, Pz_
	
	mlab.figure(size=(500,500),fgcolor=(0,0,0),bgcolor=(1,1,1))
	mlab.pipeline.volume(mlab.pipeline.scalar_field(E), vmin=0, vmax=1)
	mlab.outline(extent=[0,200,0,200,0,200],opacity=0.2,color=(0,0,0))
	mlab.text(0.65,0.9,'t = '+'%.4f'%(-1.5-0.5*Dt + n*Dt)+' s',width=0.3)
	mlab.text(0.5,0.8,'Max. teor. = '+'%.2f'%max_val,width=0.45)

mlab.show()
