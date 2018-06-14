import numpy as np
import h5py as h5
from mayavi import mlab
from F25 import *
import sys

SEEDS = 8
n=int(sys.argv[3])
ARG1 = sys.argv[1]
ARG2 = sys.argv[2]
try:
	ARG3 = sys.argv[4]
except:
	ARG3 = None

def ind(i,j,k):
	global N
	return (N+1)*(N+1)*i+(N+1)*j+k


if ARG2 == 'sim':
	###### INPUT DATA ######
	
	fin = h5.File('results_25.h5','r')
	N = fin.attrs['N'][0]
	Dt = fin.attrs['Dt'][0]
	
	I = np.array([[[i for k in range(N)] for j in range(N)] for i in range(N)])+1
	J = np.array([[[j for k in range(N)] for j in range(N)] for i in range(N)])+1
	K = np.array([[[k for k in range(N)] for j in range(N)] for i in range(N)])+1
	
	if ARG3 != None:
		mlab.figure(size=(700,700),fgcolor=(0,0,0),bgcolor=(1,1,1))
	
	else:
		mlab.figure(size=(1366,1366),fgcolor=(0,0,0),bgcolor=(1,1,1))
	
	if ARG1 == 'E':
		Ex = fin[str(n)]['Ex'][:]
		Ey = fin[str(n)]['Ey'][:]
		Ez = fin[str(n)]['Ez'][:]
		
		Ex = 0.5  *( Ex[ind(I-1,J-1,K-1)] + Ex[ind( I ,J-1,K-1)] )
		Ey = 0.5  *( Ey[ind(I-1,J-1,K-1)] + Ey[ind(I-1, J ,K-1)] )
		Ez = 0.5  *( Ez[ind(I-1,J-1,K-1)] + Ez[ind(I-1,J-1, K )] )
		
		field = mlab.pipeline.extract_vector_norm(mlab.pipeline.vector_field(Ex,Ey,Ez)); del Ex, Ey, Ez;
		mlab.text(0.65,0.9,'t = '+"%.4f"%(-1.5-0.5*Dt + n*Dt)+' s',width=0.3)
		color = (30./255,144./255,1)
		
		mlab.outline(extent=[0,200,0,200,0,200],opacity=0.2,color=(0,0,0))
		mlab.orientation_axes()
		
		if n == 60:
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [100.5,100.5,100.75]
			line.seed.widget.radius = 10
			line.seed.widget.theta_resolution=5
			line.seed.widget.phi_resolution=5
			
			#line.seed.widget.enabled=False
		
		elif n == 1:
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [100,100,100]
			line.seed.widget.radius = 10
			line.seed.widget.theta_resolution=5
			line.seed.widget.phi_resolution=5
			
			#line.seed.widget.enabled=False
		
		elif n == 120:
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [100,100,100]
			line.seed.widget.radius = 10
			line.seed.widget.theta_resolution=5
			line.seed.widget.phi_resolution=5
			
			#line.seed.widget.enabled=False
	
	
	elif ARG1 == 'H':
		Hx = fin[str(n)]['Hx'][:]
		Hy = fin[str(n)]['Hy'][:]
		Hz = fin[str(n)]['Hz'][:]
		
		Hx = 0.5  *( Hx[ind(I,J,K)] + Hx[ind(I-1, J , K )] )
		Hy = 0.5  *( Hy[ind(I,J,K)] + Hy[ind( I ,J-1, K )] )
		Hz = 0.5  *( Hz[ind(I,J,K)] + Hz[ind( I , J ,K-1)] )
		
		field = mlab.pipeline.extract_vector_norm(mlab.pipeline.vector_field(Hx,Hy,Hz)); del Hx, Hy, Hz;
		mlab.text(0.65,0.9,'t = '+"%.4f"%(-1.5 + n*Dt)+' s',width=0.3)
		color = (1,165./255,0)
		
		mlab.outline(extent=[0,200,0,200,0,200],opacity=0.2,color=(0,0,0))
		mlab.orientation_axes()
		
		if n == 60:
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [100,100,100]
			line.seed.widget.radius = 10
			line.seed.widget.theta_resolution=5
			line.seed.widget.phi_resolution=5
			
			#line.seed.widget.enabled=False
		
		elif n == 1:
			
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [100,100,100]
			line.seed.widget.radius = 10
			line.seed.widget.theta_resolution=5
			line.seed.widget.phi_resolution=5
			
			#line.seed.widget.enabled=False
		
		elif n == 120:
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [100,100,100]
			line.seed.widget.radius = 10
			line.seed.widget.theta_resolution=5
			line.seed.widget.phi_resolution=5
			
			#line.seed.widget.enabled=False
	
	
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
		mlab.text(0.65,0.9,'t = '+"%.4f"%(-1.5-0.5*Dt + n*Dt)+' s',width=0.3)
		color = (0,1,0)
		
		mlab.outline(extent=[0,200,0,200,0,200],opacity=0.2,color=(0,0,0))
		mlab.orientation_axes()
		
		mlab.pipeline.iso_surface(module, contours=[P.max()-0.5*P.ptp()], opacity=0.5,reset_zoom=False,color=(1,0,0))
	
		if n == 60:
			axe = mlab.pipeline.streamline(field,seedtype='point',seed_visible=True,integration_direction='both',color=(0,1./2.55,0))
			axe.seed.widget.position=[100,100,100]
			axe.streamline_type='tube'
			axe.tube_filter.radius = 3
			axe.stream_tracer.maximum_propagation = 200
			axe.seed.widget.enabled=False
			
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [100,100,100]
			line.seed.widget.radius = 40
			line.seed.widget.theta_resolution=20
			line.seed.widget.phi_resolution=2
			
			line.seed.widget.enabled=False
		
		elif n == 1:
			
			axe = mlab.pipeline.streamline(field,seedtype='point',seed_visible=True,integration_direction='both',color=(0,1./2.55,0))
			axe.seed.widget.position=[100,100,100]
			axe.streamline_type='tube'
			axe.tube_filter.radius = 3
			axe.stream_tracer.maximum_propagation = 200
			axe.seed.widget.enabled=False
			
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [100,100,70.25]
			line.seed.widget.radius = 40
			line.seed.widget.theta_resolution=20
			line.seed.widget.phi_resolution=2
			
			line.seed.widget.enabled=False
		
		elif n == 120:
			axe = mlab.pipeline.streamline(field,seedtype='point',seed_visible=True,integration_direction='both',color=(0,1./2.55,0))
			axe.seed.widget.position=[100,100,100]
			axe.streamline_type='tube'
			axe.tube_filter.radius = 3
			axe.stream_tracer.maximum_propagation = 200
			axe.seed.widget.enabled=False
			
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [100,100,129.93]
			line.seed.widget.radius = 40
			line.seed.widget.theta_resolution=20
			line.seed.widget.phi_resolution=2
			
			line.seed.widget.enabled=False
	
	
	elif ARG1 == 'EHP':
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
		
		fieldP = mlab.pipeline.extract_vector_norm(mlab.pipeline.vector_field(Px,Py,Pz)); del Px, Py, Pz;
		fieldE = mlab.pipeline.extract_vector_norm(mlab.pipeline.vector_field(Ex,Ey,Ez)); del Ex, Ey, Ez;
		fieldH = mlab.pipeline.extract_vector_norm(mlab.pipeline.vector_field(Hx,Hy,Hz)); del Hx, Hy, Hz;
		mlab.text(0.65,0.9,'t = '+"%.4f"%(-1.5-0.5*Dt + n*Dt)+' s',width=0.3)
		colorP = (0,1,0)
		colorE = (30./255,144./255,1)
		colorH = (1,165./255,0)
		
		mlab.outline(extent=[0,200,0,200,0,200],opacity=0.2,color=(0,0,0))
		mlab.orientation_axes()
		
		if n == 60:
			axeP = mlab.pipeline.streamline(fieldP,seedtype='point',seed_visible=True,integration_direction='both',color=(0,1./2.55,0))
			axeP.seed.widget.position=[100,100,100]
			axeP.streamline_type='tube'
			axeP.tube_filter.radius = 3
			axeP.stream_tracer.maximum_propagation = 200
			axeP.seed.widget.enabled=False
			
			axeH = mlab.pipeline.streamline(fieldH,seedtype='point',seed_visible=True,integration_direction='both',color=(1,69./255,0))
			axeH.seed.widget.position=[100,100,100.25]
			axeH.streamline_type='tube'
			axeH.tube_filter.radius = 3
			axeH.stream_tracer.maximum_propagation = 200
			axeH.seed.widget.enabled=False
			
			axeE = mlab.pipeline.streamline(fieldE,seedtype='point',seed_visible=True,integration_direction='both',color=(0,0,139./255))
			axeE.seed.widget.position=[100,100,100.25]
			axeE.streamline_type='tube'
			axeE.tube_filter.radius = 3
			axeE.stream_tracer.maximum_propagation = 200
			axeE.seed.widget.enabled=False
			
			lineP = mlab.pipeline.streamline(fieldP, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorP)
			lineP.streamline_type = 'tube'
			lineP.tube_filter.radius = 2
			lineP.stream_tracer.maximum_propagation = 200
			lineP.seed.widget.center = [100,100,100]
			lineP.seed.widget.radius = 40
			lineP.seed.widget.theta_resolution=20
			lineP.seed.widget.phi_resolution=2
			
			lineH = []
			for i in range(SEEDS):
				lineH += [ mlab.pipeline.streamline(fieldH, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorH) ]
				lineH[i].streamline_type = 'tube'
				lineH[i].tube_filter.radius = 2
				lineH[i].stream_tracer.maximum_propagation = 200
				lineH[i].seed.widget.center = [100+30.5*np.cos(1.*i/SEEDS*2*np.pi),100,100+30.5*np.sin(1.*i/SEEDS*2*np.pi)]
				lineH[i].seed.widget.radius = 10
				lineH[i].seed.widget.theta_resolution=5
				lineH[i].seed.widget.phi_resolution=5
			
			lineE = []
			for i in range(SEEDS):
				lineE += [ mlab.pipeline.streamline(fieldE, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorE) ]
				lineE[i].streamline_type = 'tube'
				lineE[i].tube_filter.radius = 2
				lineE[i].stream_tracer.maximum_propagation = 200
				lineE[i].seed.widget.center = [100,100+30.5*np.cos(1.*i/SEEDS*2*np.pi),100+30.5*np.sin(1.*i/SEEDS*2*np.pi)]
				lineE[i].seed.widget.radius = 10
				lineE[i].seed.widget.theta_resolution=5
				lineE[i].seed.widget.phi_resolution=5
			
			lineP.seed.widget.enabled=False
			for i in lineH:
				i.seed.widget.enabled=False
			for i in lineE:
				i.seed.widget.enabled=False
		
		elif n == 1:
			axeP = mlab.pipeline.streamline(fieldP,seedtype='point',seed_visible=True,integration_direction='both',color=(0,1./2.55,0))
			axeP.seed.widget.position=[100,100,100]
			axeP.streamline_type='tube'
			axeP.tube_filter.radius = 3
			axeP.stream_tracer.maximum_propagation = 200
			axeP.seed.widget.enabled=False
			
			axeH = mlab.pipeline.streamline(fieldH,seedtype='point',seed_visible=True,integration_direction='both',color=(1,69./255,0))
			axeH.seed.widget.position=[100,100,129.75]
			axeH.streamline_type='tube'
			axeH.tube_filter.radius = 3
			axeH.stream_tracer.maximum_propagation = 200
			axeH.seed.widget.enabled=False
			
			axeE = mlab.pipeline.streamline(fieldE,seedtype='point',seed_visible=True,integration_direction='both',color=(0,0,139./255))
			axeE.seed.widget.position=[100,100,129.75]
			axeE.streamline_type='tube'
			axeE.tube_filter.radius = 3
			axeE.stream_tracer.maximum_propagation = 200
			axeE.seed.widget.enabled=False
			
			lineP = mlab.pipeline.streamline(fieldP, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorP)
			lineP.streamline_type = 'tube'
			lineP.tube_filter.radius = 2
			lineP.stream_tracer.maximum_propagation = 200
			lineP.seed.widget.center = [100,100,70.13]
			lineP.seed.widget.radius = 40
			lineP.seed.widget.theta_resolution=20
			lineP.seed.widget.phi_resolution=2
			
			lineH = mlab.pipeline.streamline(fieldH, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorH)
			lineH.streamline_type = 'tube'
			lineH.tube_filter.radius = 2
			lineH.stream_tracer.maximum_propagation = 200
			lineH.seed.widget.center = [100,100,100]
			lineH.seed.widget.radius = 20
			lineH.seed.widget.theta_resolution=15
			lineH.seed.widget.phi_resolution=15
			
			lineE = mlab.pipeline.streamline(fieldE, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorE)
			lineE.streamline_type = 'tube'
			lineE.tube_filter.radius = 2
			lineE.stream_tracer.maximum_propagation = 200
			lineE.seed.widget.center = [100,100,100]
			lineE.seed.widget.radius = 20
			lineE.seed.widget.theta_resolution=15
			lineE.seed.widget.phi_resolution=15
			
			lineP.seed.widget.enabled=False
			lineH.seed.widget.enabled=False
			lineE.seed.widget.enabled=False
		
		elif n == 120:
			axeP = mlab.pipeline.streamline(fieldP,seedtype='point',seed_visible=True,integration_direction='both',color=(0,1./2.55,0))
			axeP.seed.widget.position=[100,100,100]
			axeP.streamline_type='tube'
			axeP.tube_filter.radius = 3
			axeP.stream_tracer.maximum_propagation = 200
			axeP.seed.widget.enabled=False
			
			axeH = mlab.pipeline.streamline(fieldH,seedtype='point',seed_visible=True,integration_direction='both',color=(1,69./255,0))
			axeH.seed.widget.position=[100,100,70.25]
			axeH.streamline_type='tube'
			axeH.tube_filter.radius = 3
			axeH.stream_tracer.maximum_propagation = 200
			axeH.seed.widget.enabled=False
			
			axeE = mlab.pipeline.streamline(fieldE,seedtype='point',seed_visible=True,integration_direction='both',color=(0,0,139./255))
			axeE.seed.widget.position=[100,100,70.25]
			axeE.streamline_type='tube'
			axeE.tube_filter.radius = 3
			axeE.stream_tracer.maximum_propagation = 200
			axeE.seed.widget.enabled=False
			
			lineP = mlab.pipeline.streamline(fieldP, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorP)
			lineP.streamline_type = 'tube'
			lineP.tube_filter.radius = 2
			lineP.stream_tracer.maximum_propagation = 200
			lineP.seed.widget.center = [100,100,129.93]
			lineP.seed.widget.radius = 40
			lineP.seed.widget.theta_resolution=20
			lineP.seed.widget.phi_resolution=2
			
			lineH = mlab.pipeline.streamline(fieldH, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorH)
			lineH.streamline_type = 'tube'
			lineH.tube_filter.radius = 2
			lineH.stream_tracer.maximum_propagation = 200
			lineH.seed.widget.center = [100,100,100]
			lineH.seed.widget.radius = 20
			lineH.seed.widget.theta_resolution=15
			lineH.seed.widget.phi_resolution=15
			
			lineE = mlab.pipeline.streamline(fieldE, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorE)
			lineE.streamline_type = 'tube'
			lineE.tube_filter.radius = 2
			lineE.stream_tracer.maximum_propagation = 200
			lineE.seed.widget.center = [100,100,100]
			lineE.seed.widget.radius = 20
			lineE.seed.widget.theta_resolution=15
			lineE.seed.widget.phi_resolution=15
			
			lineP.seed.widget.enabled=False
			lineH.seed.widget.enabled=False
			lineE.seed.widget.enabled=False
	
	'''
	if ARG3 == None:
		print 'caca'
		mlab.view(distance=700)
		if n == 1:
			mlab.savefig('../hopfions_memoria/media/23_flow_t-15_'+ARG1+'_'+ARG2+'.png',size=(1366,1366))
		elif n == 60:
			mlab.savefig('../hopfions_memoria/media/23_flow_t0_'+ARG1+'_'+ARG2+'.png',size=(1366,1366))
		elif n == 120:
			mlab.savefig('../hopfions_memoria/media/23_flow_t15_'+ARG1+'_'+ARG2+'.png',size=(1366,1366))
	
	else:
		mlab.view(azimuth=0,elevation=0,distance=600)
		if n == 1:
			mlab.savefig('../hopfions_memoria/media/23_flow_t-15_'+ARG1+'_'+ARG2+'_above.png',size=(700,700))
		elif n == 60:
			mlab.savefig('../hopfions_memoria/media/23_flow_t0_'+ARG1+'_'+ARG2+'_above.png',size=(700,700))
		elif n == 120:
			mlab.savefig('../hopfions_memoria/media/23_flow_t15_'+ARG1+'_'+ARG2+'_above.png',size=(700,700))
	'''
	mlab.show()

###################################################################### ^ SIMULADO ^


#                        ~Territorio de nadie~						 #


###################################################################### v TEORICO  v

elif ARG2 == 'teor':
	###### INPUT DATA ######
	
	N=200
	D = 10./N
	Dt = 0.5*D
	
	I = D * (np.array([[[i for k in range(N+1)] for j in range(N+1)] for i in range(N+1)]) - 0.5*N)
	J = D * (np.array([[[j for k in range(N+1)] for j in range(N+1)] for i in range(N+1)]) - 0.5*N)
	K = D * (np.array([[[k for k in range(N+1)] for j in range(N+1)] for i in range(N+1)]) - 0.5*N)
	
	if ARG3 != None:
		mlab.figure(size=(700,700),fgcolor=(0,0,0),bgcolor=(1,1,1))
	
	else:
		mlab.figure(size=(1366,1366),fgcolor=(0,0,0),bgcolor=(1,1,1))
	
	if ARG1 == 'E':
		Ex = ( Fx( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		Ey = ( Fy( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		Ez = ( Fz( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		
		field = mlab.pipeline.extract_vector_norm(mlab.pipeline.vector_field(Ex,Ey,Ez)); del Ex, Ey, Ez;
		mlab.text(0.65,0.9,'t = '+"%.4f"%(-1.5-0.5*Dt + n*Dt)+' s',width=0.3)
		isosurface=False
		color = (30./255,144./255,1)
		
		mlab.outline(extent=[0,200,0,200,0,200],opacity=0.2,color=(0,0,0))
		mlab.orientation_axes()
		
		if n == 60:
			
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [100.5,100.5,100.75]
			line.seed.widget.radius = 10
			line.seed.widget.theta_resolution=5
			line.seed.widget.phi_resolution=5
			
			#line.seed.widget.enabled=False
		
		elif n == 1:
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,seed_resolution=3,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [100.5,100.5,100.5]
			line.seed.widget.radius = 10
			line.seed.widget.theta_resolution=5
			line.seed.widget.phi_resolution=5
			
			#line.seed.widget.enabled=False
			
		elif n == 120:
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [101,101,101]
			line.seed.widget.radius = 10
			line.seed.widget.theta_resolution=5
			line.seed.widget.phi_resolution=5
			
			#line.seed.widget.enabled=False
	
	
	
	elif ARG1 == 'H':
		Hx = ( Fx( I , J , K , -1.5+n*Dt ) ).imag
		Hy = ( Fy( I , J , K , -1.5+n*Dt ) ).imag
		Hz = ( Fz( I , J , K , -1.5+n*Dt ) ).imag
		
		field = mlab.pipeline.extract_vector_norm(mlab.pipeline.vector_field(Hx,Hy,Hz)); del Hx, Hy, Hz;
		mlab.text(0.65,0.9,'t = '+"%.4f"%(-1.5 + n*Dt)+' s',width=0.3)
		isosurface=False
		color = (1,165./255,0)
		
		mlab.outline(extent=[0,200,0,200,0,200],opacity=0.2,color=(0,0,0))
		mlab.orientation_axes()
		
		if n == 60:
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [101,101,101]
			line.seed.widget.radius = 10
			line.seed.widget.theta_resolution=5
			line.seed.widget.phi_resolution=5
			
			#line.seed.widget.enabled=False
		
		elif n == 1:
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [101,101,101]
			line.seed.widget.radius = 10
			line.seed.widget.theta_resolution=5
			line.seed.widget.phi_resolution=5
			
			#line.seed.widget.enabled=False
		
		elif n == 120:
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [101,101,101]
			line.seed.widget.radius = 10
			line.seed.widget.theta_resolution=5
			line.seed.widget.phi_resolution=5
			
			#line.seed.widget.enabled=False
	
	
	
	elif ARG1 == 'P':
		Ex = ( Fx( I , J , K , -1.5-0.5*Dt+n*Dt ) ).real
		Ey = ( Fy( I , J , K , -1.5-0.5*Dt+n*Dt ) ).real
		Ez = ( Fz( I , J , K , -1.5-0.5*Dt+n*Dt ) ).real
		
		Hx = ( Fx( I , J , K , -1.5-0.5*Dt+n*Dt ) ).imag
		Hy = ( Fy( I , J , K , -1.5-0.5*Dt+n*Dt ) ).imag
		Hz = ( Fz( I , J , K , -1.5-0.5*Dt+n*Dt ) ).imag
		
		Px = Ey*Hz-Ez*Hy
		Py = Ez*Hx-Ex*Hz
		Pz = Ex*Hy-Ey*Hx
		
		P = np.sqrt(Px*Px+Py*Py+Pz*Pz)
		module = mlab.pipeline.scalar_field(P)
		field = mlab.pipeline.extract_vector_norm(mlab.pipeline.vector_field(Px,Py,Pz)); del Px, Py, Pz;
		mlab.text(0.65,0.9,'t = '+"%.4f"%(-1.5-0.5*Dt + n*Dt)+' s',width=0.3)
		color = (0,1,0)
		
		mlab.outline(extent=[0,200,0,200,0,200],opacity=0.2,color=(0,0,0))
		mlab.orientation_axes()
		
		mlab.pipeline.iso_surface(module, contours=[P.max()-0.5*P.ptp()], opacity=0.5,reset_zoom=False,color=(1,0,0))
		
		if n == 60:
			axe = mlab.pipeline.streamline(field,seedtype='point',seed_visible=True,integration_direction='both',color=(0,1./2.55,0))
			axe.seed.widget.position=[101,101,101]
			axe.streamline_type='tube'
			axe.tube_filter.radius = 3
			axe.stream_tracer.maximum_propagation = 200
			axe.seed.widget.enabled=False
			
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [101,101,101]
			line.seed.widget.radius = 40
			line.seed.widget.theta_resolution=20
			line.seed.widget.phi_resolution=2
			
			line.seed.widget.enabled=False
		
		elif n == 1:
			axe = mlab.pipeline.streamline(field,seedtype='point',seed_visible=True,integration_direction='both',color=(0,1./2.55,0))
			axe.seed.widget.position=[101,101,101]
			axe.streamline_type='tube'
			axe.tube_filter.radius = 3
			axe.stream_tracer.maximum_propagation = 200
			axe.seed.widget.enabled=False
			
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [101,101,129.93]
			line.seed.widget.radius = 40
			line.seed.widget.theta_resolution=20
			line.seed.widget.phi_resolution=2
			
			line.seed.widget.enabled=False
		
		elif n == 120:
			axe = mlab.pipeline.streamline(field,seedtype='point',seed_visible=True,integration_direction='both',color=(0,1./2.55,0))
			axe.seed.widget.position=[101,101,101]
			axe.streamline_type='tube'
			axe.tube_filter.radius = 3
			axe.stream_tracer.maximum_propagation = 200
			axe.seed.widget.enabled=False
			
			line = mlab.pipeline.streamline(field, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=color)
			line.streamline_type = 'tube'
			line.tube_filter.radius = 2
			line.stream_tracer.maximum_propagation = 500
			line.seed.widget.center = [101,101,129.93]
			line.seed.widget.radius = 40
			line.seed.widget.theta_resolution=20
			line.seed.widget.phi_resolution=2
			
			line.seed.widget.enabled=False
	
	elif ARG1 == 'EHP':
		Ex = ( Fx( I , J , K , -1.5-0.5*Dt+n*Dt ) ).real
		Ey = ( Fy( I , J , K , -1.5-0.5*Dt+n*Dt ) ).real
		Ez = ( Fz( I , J , K , -1.5-0.5*Dt+n*Dt ) ).real
		
		Hx = ( Fx( I , J , K , -1.5-0.5*Dt+n*Dt ) ).imag
		Hy = ( Fy( I , J , K , -1.5-0.5*Dt+n*Dt ) ).imag
		Hz = ( Fz( I , J , K , -1.5-0.5*Dt+n*Dt ) ).imag
		
		Px = Ey*Hz-Ez*Hy
		Py = Ez*Hx-Ex*Hz
		Pz = Ex*Hy-Ey*Hx
		
		fieldP = mlab.pipeline.extract_vector_norm(mlab.pipeline.vector_field(Px,Py,Pz)); del Px, Py, Pz;
		fieldE = mlab.pipeline.extract_vector_norm(mlab.pipeline.vector_field(Ex,Ey,Ez)); del Ex, Ey, Ez;
		fieldH = mlab.pipeline.extract_vector_norm(mlab.pipeline.vector_field(Hx,Hy,Hz)); del Hx, Hy, Hz;
		mlab.text(0.65,0.9,'t = '+"%.4f"%(-1.5-0.5*Dt + n*Dt)+' s',width=0.3)
		colorP = (0,1,0)
		colorE = (30./255,144./255,1)
		colorH = (1,165./255,0)
		
		mlab.outline(extent=[0,200,0,200,0,200],opacity=0.2,color=(0,0,0))
		mlab.orientation_axes()
		
		if n == 60:
			axeP = mlab.pipeline.streamline(fieldP,seedtype='point',seed_visible=True,integration_direction='both',color=(0,1./2.55,0))
			axeP.seed.widget.position=[101,101,101]
			axeP.streamline_type='tube'
			axeP.tube_filter.radius = 3
			axeP.stream_tracer.maximum_propagation = 200
			axeP.seed.widget.enabled=False
			
			axeH = mlab.pipeline.streamline(fieldH,seedtype='point',seed_visible=True,integration_direction='both',color=(1,69./255,0))
			axeH.seed.widget.position=[101,101,101.25]
			axeH.streamline_type='tube'
			axeH.tube_filter.radius = 3
			axeH.stream_tracer.maximum_propagation = 200
			axeH.seed.widget.enabled=False
			
			axeE = mlab.pipeline.streamline(fieldE,seedtype='point',seed_visible=True,integration_direction='both',color=(0,0,139./255))
			axeE.seed.widget.position=[101,101,101.25]
			axeE.streamline_type='tube'
			axeE.tube_filter.radius = 3
			axeE.stream_tracer.maximum_propagation = 200
			axeE.seed.widget.enabled=False
			
			lineP = mlab.pipeline.streamline(fieldP, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorP)
			lineP.streamline_type = 'tube'
			lineP.tube_filter.radius = 2
			lineP.stream_tracer.maximum_propagation = 200
			lineP.seed.widget.center = [101,101,101]
			lineP.seed.widget.radius = 40
			lineP.seed.widget.theta_resolution=20
			lineP.seed.widget.phi_resolution=2
			
			lineH = []
			for i in range(SEEDS):
				lineH += [ mlab.pipeline.streamline(fieldH, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorH) ]
				lineH[i].streamline_type = 'tube'
				lineH[i].tube_filter.radius = 2
				lineH[i].stream_tracer.maximum_propagation = 200
				lineH[i].seed.widget.center = [101+30.5*np.cos(1.*i/SEEDS*2*np.pi),101,101+30.5*np.sin(1.*i/SEEDS*2*np.pi)]
				lineH[i].seed.widget.radius = 10
				lineH[i].seed.widget.theta_resolution=5
				lineH[i].seed.widget.phi_resolution=5
			
			lineE = []
			for i in range(SEEDS):
				lineE += [ mlab.pipeline.streamline(fieldE, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorE) ]
				lineE[i].streamline_type = 'tube'
				lineE[i].tube_filter.radius = 2
				lineE[i].stream_tracer.maximum_propagation = 200
				lineE[i].seed.widget.center = [101,101+30.5*np.cos(1.*i/SEEDS*2*np.pi),101+30.5*np.sin(1.*i/SEEDS*2*np.pi)]
				lineE[i].seed.widget.radius = 10
				lineE[i].seed.widget.theta_resolution=5
				lineE[i].seed.widget.phi_resolution=5
			
			lineP.seed.widget.enabled=False
			for i in lineH:
				i.seed.widget.enabled=False
			for i in lineE:
				i.seed.widget.enabled=False
		
		elif n == 1:
			axeP = mlab.pipeline.streamline(fieldP,seedtype='point',seed_visible=True,integration_direction='both',color=(0,1./2.55,0))
			axeP.seed.widget.position=[101,101,101]
			axeP.streamline_type='tube'
			axeP.tube_filter.radius = 3
			axeP.stream_tracer.maximum_propagation = 200
			axeP.seed.widget.enabled=False
			
			axeH = mlab.pipeline.streamline(fieldH,seedtype='point',seed_visible=True,integration_direction='both',color=(1,69./255,0))
			axeH.seed.widget.position=[101,101,130.75]
			axeH.streamline_type='tube'
			axeH.tube_filter.radius = 3
			axeH.stream_tracer.maximum_propagation = 200
			axeH.seed.widget.enabled=False
			
			axeE = mlab.pipeline.streamline(fieldE,seedtype='point',seed_visible=True,integration_direction='both',color=(0,0,139./255))
			axeE.seed.widget.position=[101,101,130.75]
			axeE.streamline_type='tube'
			axeE.tube_filter.radius = 3
			axeE.stream_tracer.maximum_propagation = 200
			axeE.seed.widget.enabled=False
			
			lineP = mlab.pipeline.streamline(fieldP, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorP)
			lineP.streamline_type = 'tube'
			lineP.tube_filter.radius = 2
			lineP.stream_tracer.maximum_propagation = 200
			lineP.seed.widget.center = [101,101,70.25]
			lineP.seed.widget.radius = 40
			lineP.seed.widget.theta_resolution=20
			lineP.seed.widget.phi_resolution=2
			
			lineH = mlab.pipeline.streamline(fieldH, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorH)
			lineH.streamline_type = 'tube'
			lineH.tube_filter.radius = 2
			lineH.stream_tracer.maximum_propagation = 200
			lineH.seed.widget.center = [101,101,101]
			lineH.seed.widget.radius = 20
			lineH.seed.widget.theta_resolution=15
			lineH.seed.widget.phi_resolution=15
			
			lineE = mlab.pipeline.streamline(fieldE, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorE)
			lineE.streamline_type = 'tube'
			lineE.tube_filter.radius = 2
			lineE.stream_tracer.maximum_propagation = 200
			lineE.seed.widget.center = [101,101,101]
			lineE.seed.widget.radius = 20
			lineE.seed.widget.theta_resolution=15
			lineE.seed.widget.phi_resolution=15
			
			lineP.seed.widget.enabled=False
			lineH.seed.widget.enabled=False
			lineE.seed.widget.enabled=False
		
		elif n == 120:
			axeP = mlab.pipeline.streamline(fieldP,seedtype='point',seed_visible=True,integration_direction='both',color=(0,1./2.55,0))
			axeP.seed.widget.position=[101,101,101]
			axeP.streamline_type='tube'
			axeP.tube_filter.radius = 3
			axeP.stream_tracer.maximum_propagation = 200
			axeP.seed.widget.enabled=False
			
			axeH = mlab.pipeline.streamline(fieldH,seedtype='point',seed_visible=True,integration_direction='both',color=(1,69./255,0))
			axeH.seed.widget.position=[101,101,71.25]
			axeH.streamline_type='tube'
			axeH.tube_filter.radius = 3
			axeH.stream_tracer.maximum_propagation = 200
			axeH.seed.widget.enabled=False
			
			axeE = mlab.pipeline.streamline(fieldE,seedtype='point',seed_visible=True,integration_direction='both',color=(0,0,139./255))
			axeE.seed.widget.position=[101,101,71.25]
			axeE.streamline_type='tube'
			axeE.tube_filter.radius = 3
			axeE.stream_tracer.maximum_propagation = 200
			axeE.seed.widget.enabled=False
			
			lineP = mlab.pipeline.streamline(fieldP, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorP)
			lineP.streamline_type = 'tube'
			lineP.tube_filter.radius = 2
			lineP.stream_tracer.maximum_propagation = 200
			lineP.seed.widget.center = [101,101,129.93]
			lineP.seed.widget.radius = 40
			lineP.seed.widget.theta_resolution=20
			lineP.seed.widget.phi_resolution=2
			
			lineH = mlab.pipeline.streamline(fieldH, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorH)
			lineH.streamline_type = 'tube'
			lineH.tube_filter.radius = 2
			lineH.stream_tracer.maximum_propagation = 200
			lineH.seed.widget.center = [101,101,101]
			lineH.seed.widget.radius = 20
			lineH.seed.widget.theta_resolution=15
			lineH.seed.widget.phi_resolution=15
			
			lineE = mlab.pipeline.streamline(fieldE, seedtype='sphere',seed_visible=True,integration_direction='both',vmin=0,vmax=1,color=colorE)
			lineE.streamline_type = 'tube'
			lineE.tube_filter.radius = 2
			lineE.stream_tracer.maximum_propagation = 200
			lineE.seed.widget.center = [101,101,101]
			lineE.seed.widget.radius = 20
			lineE.seed.widget.theta_resolution=15
			lineE.seed.widget.phi_resolution=15
			
			lineP.seed.widget.enabled=False
			lineH.seed.widget.enabled=False
			lineE.seed.widget.enabled=False
	
	'''
	if ARG3 ==None:
		mlab.view(distance=700)
		if n == 1:
			mlab.savefig('../hopfions_memoria/media/23_flow_t-15_'+ARG1+'_'+ARG2+'.png',size=(1366,1366))
		elif n == 60:
			mlab.savefig('../hopfions_memoria/media/23_flow_t0_'+ARG1+'_'+ARG2+'.png',size=(1366,1366))
		elif n == 120:
			mlab.savefig('../hopfions_memoria/media/23_flow_t15_'+ARG1+'_'+ARG2+'.png',size=(1366,1366))
	
	else:
		mlab.view(azimuth=0,elevation=0,distance=600)
		if n == 1:
			mlab.savefig('../hopfions_memoria/media/23_flow_t-15_'+ARG1+'_'+ARG2+'_above.png',size=(700,700))
		elif n == 60:
			mlab.savefig('../hopfions_memoria/media/23_flow_t0_'+ARG1+'_'+ARG2+'_above.png',size=(700,700))
		elif n == 120:
			mlab.savefig('../hopfions_memoria/media/23_flow_t15_'+ARG1+'_'+ARG2+'_above.png',size=(700,700))
	'''
	mlab.show()

else:
	print 'ARGUMENT ERROR'