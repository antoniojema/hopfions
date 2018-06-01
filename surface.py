import numpy as np
import h5py as h5
from mayavi import mlab
from mayavi import tools

def ind(i,j,k):
	global N
	return (N+1)*(N+1)*i+(N+1)*j+k

fin = h5.File('simulation_results.h5','r')
N = fin.attrs['N'][0]
iterations = fin.attrs['iterations'][0]
Dt = fin.attrs['Dt'][0]

I = np.array([[[i for k in range(N)] for j in range(N)] for i in range(N)])
J = np.array([[[j for k in range(N)] for j in range(N)] for i in range(N)])
K = np.array([[[k for k in range(N)] for j in range(N)] for i in range(N)])

### CODIGO CHUNGO ###
print 0

Ex = fin['0']['Ex'][:]
Ey = fin['0']['Ey'][:]
Ez = fin['0']['Ez'][:]

Ex = 0.5 * ( Ex[ind(I,J,K)] + Ex[ind(I+1, J , K )] )	
Ey = 0.5 * ( Ey[ind(I,J,K)] + Ey[ind( I ,J+1, K )] )
Ez = 0.5 * ( Ez[ind(I,J,K)] + Ez[ind( I , J ,K+1)] )
#P = Ex*Ex + Ey*Ey + Ez*Ez

#P_ = mlab.pipeline.scalar_field(P)
E = mlab.pipeline.vector_field(Ex,Ey,Ez)
E_norm = mlab.pipeline.extract_vector_norm(E)

mlab.outline(extent=[0,len(Ex),0,len(Ex),0,len(Ex)])

#surf1 = mlab.pipeline.iso_surface(src, contours=[P.max()-0.2*P.ptp()], opacity=0.5,reset_zoom=False)

#vol = mlab.pipeline.volume(src,vmin=0,vmax=0.8)

line = tools.pipeline.streamline(E_norm,seed_scale=1,seed_resolution=2,integration_direction='forward',seed_visible=True, seedtype='point', extent=[0.1*len(Ex),0.9*len(Ex),0.1*len(Ex),0.9*len(Ex),0.1*len(Ex),0.9*len(Ex)])
line.stream_tracer.maximum_propagation = 1000
'''
text = mlab.text(0.1,0.8,'t = '+"%.3f"%(-1.5-0.5*Dt + 0*Dt)+' s',width=0.4)
mlab.savefig('isosurface/P02'+'%.3i'%0+'.png')
for i in np.arange(81,iterations,1):
	print i
	text.remove()
	text = mlab.text(0.1,0.8,'t = '+"%.3f"%(-1.5-0.5*Dt + i*Dt)+' s',width=0.4)
	
	Ex = fin[str(i)]['Ex'][:]
	Ey = fin[str(i)]['Ey'][:]
	Ez = fin[str(i)]['Ez'][:]

	Ex = 0.5 * ( Ex[ind(I,J,K)] + Ex[ind(I+1, J , K )] )	
	Ey = 0.5 * ( Ey[ind(I,J,K)] + Ey[ind( I ,J+1, K )] )
	Ez = 0.5 * ( Ez[ind(I,J,K)] + Ez[ind( I , J ,K+1)] )
	P = Ex*Ex + Ey*Ey + Ez*Ez

	src = mlab.pipeline.scalar_field(P)
	if i%2==0:
		surf1 = mlab.pipeline.iso_surface(src, contours=[P.max()-0.2*P.ptp()], opacity=0.5,reset_zoom=False)
		surf2.remove()
	else:
		surf2 = mlab.pipeline.iso_surface(src, contours=[P.max()-0.2*P.ptp()], opacity=0.5,reset_zoom=False)
		surf1.remove()
	
	#vol = mlab.pipeline.volume(src,vmin=0,vmax=0.8)
	
	#line = tools.pipeline.streamline(Ex,Ey,Ez,seed_scale=1,seed_resolution=5,integration_direction='both',seed_visible=True, seedtype='sphere', extent=[0.1*len(Ex),0.9*len(Ex),0.1*len(Ex),0.9*len(Ex),0.1*len(Ex),0.9*len(Ex)])
	
	mlab.savefig('isosurface/P02'+'%.3i'%i+'.png')
'''
mlab.show()
