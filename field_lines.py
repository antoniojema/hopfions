import numpy as np
import h5py as h5
from mayavi import mlab
from F11 import *
import sys

n=60
ARG1 = sys.argv[1]
ARG2=  sys.argv[2]

def ind(i,j,k):
	global N
	return (N+1)*(N+1)*i+(N+1)*j+k

if ARG2 == 'sim':
	###### INPUT DATA ######
	
	fin = h5.File('results_11.h5','r')
	N = fin.attrs['N'][0]
	Dt = fin.attrs['Dt'][0]
	
	I = np.array([[[i for k in range(N)] for j in range(N)] for i in range(N)])+1
	J = np.array([[[j for k in range(N)] for j in range(N)] for i in range(N)])+1
	K = np.array([[[k for k in range(N)] for j in range(N)] for i in range(N)])+1
	
	mlab.figure(size=(280,280),fgcolor=(0,0,0),bgcolor=(1,1,1))
	
	if ARG1 == 'E':
		Ex = fin[str(n)]['Ex'][:]
		Ey = fin[str(n)]['Ey'][:]
		Ez = fin[str(n)]['Ez'][:]
		
		Ex = 0.5  *( Ex[ind(I-1,J-1,K-1)] + Ex[ind( I ,J-1,K-1)] )
		Ey = 0.5  *( Ey[ind(I-1,J-1,K-1)] + Ey[ind(I-1, J ,K-1)] )
		Ez = 0.5  *( Ez[ind(I-1,J-1,K-1)] + Ez[ind(I-1,J-1, K )] )
		
		field = mlab.pipeline.extract_vector_norm(mlab.pipeline.vector_field(Ex,Ey,Ez)); del Ex, Ey, Ez;
		mlab.view(azimuth=70,elevation=70,distance=600,roll=131)
		mlab.text(0.05,0.9,'(a)',width=0.05)
		mlab.text(0.65,0.9,'t = '+"%.4f"%(-1.5-0.5*Dt + n*Dt)+' s',width=0.3)
		isosurface=False
		color = (0,0,1)
	
	elif ARG1 == 'H':
		Hx = fin[str(n)]['Hx'][:]
		Hy = fin[str(n)]['Hy'][:]
		Hz = fin[str(n)]['Hz'][:]
		
		Hx = 0.5  *( Hx[ind(I,J,K)] + Hx[ind(I-1, J , K )] )
		Hy = 0.5  *( Hy[ind(I,J,K)] + Hy[ind( I ,J-1, K )] )
		Hz = 0.5  *( Hz[ind(I,J,K)] + Hz[ind( I , J ,K-1)] )
		
		field = mlab.pipeline.extract_vector_norm(mlab.pipeline.vector_field(Hx,Hy,Hz)); del Hx, Hy, Hz;
		mlab.view(azimuth=40,elevation=25,distance=700,roll=0)
		mlab.text(0.05,0.9,'(b)',width=0.05)
		mlab.text(0.65,0.9,'t = '+"%.4f"%(-1.5 + n*Dt)+' s',width=0.3)
		isosurface=False
		color = (1,165./255,0)
        
	elif ARG1 == 'P':
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
		module = mlab.pipeline.scalar_field(P)
		field = mlab.pipeline.extract_vector_norm(mlab.pipeline.vector_field(Px,Py,Pz)); del Px, Py, Pz;
		mlab.view(azimuth=20,elevation=70,distance=700)
		mlab.text(0.65,0.9,'t = '+"%.4f"%(-1.5-0.5*Dt + n*Dt)+' s',width=0.3)
		isosurface=True
		color = (0.4,0.4,0.4)
	
	
	mlab.outline(extent=[0,200,0,200,0,200],opacity=0.2,color=(0,0,0))
	mlab.orientation_axes()
	#mlab.pipeline.image_plane_widget(module,plane_orientation='x_axes',slice_index=10)
	#mlab.pipeline.image_plane_widget(module,plane_orientation='y_axes',slice_index=10)
	#mlab.pipeline.image_plane_widget(module,plane_orientation='z_axes',slice_index=10)
	if isosurface:
		mlab.pipeline.iso_surface(module, contours=[P.max()-0.5*P.ptp()], opacity=0.5,reset_zoom=False,color=(1,0,0))
	
	line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1)
	line.stream_tracer.maximum_propagation = 1000
	
	line.seed.widget.center = [100,100,100]
	if isosurface:
		line.seed.widget.radius = 70
		line.seed.widget.theta_resolution=10
		line.seed.widget.phi_resolution=2
	else:
		line.seed.widget.theta_resolution=4
		line.seed.widget.phi_resolution=4
	
	line.seed.widget.enabled=False
	
	mlab.show()



elif ARG2 == 'teor':
	###### INPUT DATA ######
	
	N=200
	D = 10./N
	Dt = 0.5*D
	
	I = D * (np.array([[[i for k in range(N+1)] for j in range(N+1)] for i in range(N+1)]) - 0.5*N)
	J = D * (np.array([[[j for k in range(N+1)] for j in range(N+1)] for i in range(N+1)]) - 0.5*N)
	K = D * (np.array([[[k for k in range(N+1)] for j in range(N+1)] for i in range(N+1)]) - 0.5*N)
	
	mlab.figure(size=(500,500),fgcolor=(0,0,0),bgcolor=(1,1,1))
	
	if ARG1 == 'E':
		Ex = ( Fx( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		Ey = ( Fy( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		Ez = ( Fz( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		
		field = mlab.pipeline.extract_vector_norm(mlab.pipeline.vector_field(Ex,Ey,Ez)); del Ex, Ey, Ez;
		mlab.view(azimuth=70,elevation=70,distance=600,roll=131)
		mlab.text(0.05,0.9,'(a)',width=0.05)
		mlab.text(0.65,0.9,'t = '+"%.4f"%(-1.5-0.5*Dt + n*Dt)+' s',width=0.3)
		isosurface=False
		color = (0,0,1)
	
	elif ARG1 == 'H':
		Hx = ( Fx( I , J , K , -1.5+n*Dt ) ).imag
		Hy = ( Fy( I , J , K , -1.5+n*Dt ) ).imag
		Hz = ( Fz( I , J , K , -1.5+n*Dt ) ).imag
		
		field = mlab.pipeline.extract_vector_norm(mlab.pipeline.vector_field(Hx,Hy,Hz)); del Hx, Hy, Hz;
		mlab.view(azimuth=40,elevation=25,distance=700,roll=0)
		mlab.text(0.05,0.9,'(b)',width=0.05)
		mlab.text(0.65,0.9,'t = '+"%.4f"%(-1.5 + n*Dt)+' s',width=0.3)
		isosurface=False
		color = (1,165./255,0)
	
	elif ARG1 == 'P':
		Ex = ( Fx( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		Ey = ( Fy( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		Ez = ( Fz( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		
		Hx = ( Fx( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+n*Dt ) ).imag
		Hy = ( Fy( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+n*Dt ) ).imag
		Hz = ( Fz( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+n*Dt ) ).imag
		
		Px = Ey*Hz-Ez*Hy
		Py = Ez*Hx-Ex*Hz
		Pz = Ex*Hy-Ey*Hx
		del Ex, Ey, Ez, Hx, Hy, Hz
		
		P = np.sqrt(Px*Px+Py*Py+Pz*Pz)
		module = mlab.pipeline.scalar_field(P)
		field = mlab.pipeline.extract_vector_norm(mlab.pipeline.vector_field(Px,Py,Pz)); del Px, Py, Pz;
		mlab.view(azimuth=20,elevation=70,distance=700)
		mlab.text(0.65,0.9,'t = '+"%.4f"%(-1.5-0.5*Dt + n*Dt)+' s',width=0.3)
		color = (0.4,0.4,0.4)
		isosurface=True
	
	
	mlab.outline(extent=[0,200,0,200,0,200],opacity=0.2,color=(0,0,0))
	mlab.orientation_axes()
	#mlab.pipeline.image_plane_widget(module,plane_orientation='x_axes',slice_index=10)
	#mlab.pipeline.image_plane_widget(module,plane_orientation='y_axes',slice_index=10)
	#mlab.pipeline.image_plane_widget(module,plane_orientation='z_axes',slice_index=10)
	if isosurface:
		mlab.pipeline.iso_surface(module, contours=[P.max()-0.5*P.ptp()], opacity=0.5,reset_zoom=False,color=(1,0,0))
	
	line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1)
	line.stream_tracer.maximum_propagation = 1000
	
	line.seed.widget.center = [100,100,100]
	if isosurface:
		line.seed.widget.radius = 70
		line.seed.widget.theta_resolution=10
		line.seed.widget.phi_resolution=2
	else:
		line.seed.widget.theta_resolution=4
		line.seed.widget.phi_resolution=4
	
	line.seed.widget.enabled=False
	
	mlab.show()

else:
	print 'ARGUMENT ERROR'