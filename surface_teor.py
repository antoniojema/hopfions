import numpy as np
import time
from mayavi import mlab
from F23 import *

N=200
D = 10./N
Dt = 0.5*D
iterations = 200

I = D * (np.array([[[i for k in range(N)] for j in range(N)] for i in range(N)]) - 0.5*N)
J = D * (np.array([[[j for k in range(N)] for j in range(N)] for i in range(N)]) - 0.5*N)
K = D * (np.array([[[k for k in range(N)] for j in range(N)] for i in range(N)]) - 0.5*N)

### CODIGO CHUNGO ###
print 0

Ex = ( Fx( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
Ey = ( Fy( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
Ez = ( Fz( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
P = Ex*Ex + Ey*Ey + Ez*Ez

P_ = mlab.pipeline.scalar_field(P)
#E = mlab.pipeline.vector_field(Ex,Ey,Ez)
#E_norm = mlab.pipeline.extract_vector_norm(E)

#mlab.outline(extent=[0,len(Ex),0,len(Ex),0,len(Ex)])

surf1 = mlab.pipeline.iso_surface(P_, contours=[P.max()-0.2*P.ptp()], opacity=0.5,reset_zoom=False)

#vol = mlab.pipeline.volume(P_,vmin=0,vmax=0.8)

#line = mlab.flow(Ex,Ey,Ez,seed_scale=1,seed_resolution=5,integration_direction='both',seed_visible=True, seedtype='sphere', extent=[0.1*len(Ex),0.9*len(Ex),0.1*len(Ex),0.9*len(Ex),0.1*len(Ex),0.9*len(Ex)])

text = mlab.text(0.1,0.8,'t = '+"%.3f"%(-1.5-0.5*Dt + 0*Dt)+' s',width=0.4)
mlab.savefig('isosurface/P02_teor'+'%.3i'%0+'.png')

for i in np.arange(111,iterations,1):
	print i
	
	Ex = ( Fx( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+i*Dt ) ).real
	Ey = ( Fy( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+i*Dt ) ).real
	Ez = ( Fz( I+0.5*D , J+0.5*D , K+0.5*D , -1.5-0.5*Dt+i*Dt ) ).real
	
	P = Ex*Ex + Ey*Ey + Ez*Ez
	#E = mlab.pipeline.vector_field(Ex,Ey,Ez)
	#E_norm = mlab.pipeline.extract_vector_norm(E)
	
	P_ = mlab.pipeline.scalar_field(P)
	if i%2==0:
		surf1 = mlab.pipeline.iso_surface(P_, contours=[P.max()-0.2*P.ptp()],opacity=0.5,reset_zoom=False)
		surf2.remove()
	else:
		surf2 = mlab.pipeline.iso_surface(P_, contours=[P.max()-0.2*P.ptp()],opacity=0.5,reset_zoom=False)
		surf1.remove()
	
	#mlab.pipeline.volume(src,vmin=0,vmax=0.8)
	
	#mlab.flow(Ex,Ey,Ez,seed_scale=1,seed_resolution=5,integration_direction='both',seed_visible=True, seedtype='sphere', extent=[0.1*len(Ex),0.9*len(Ex),0.1*len(Ex),0.9*len(Ex),0.1*len(Ex),0.9*len(Ex)])
	
	text.remove()
	text = mlab.text(0.1,0.8,'t = '+"%.3f"%(-1.5-0.5*Dt + i*Dt)+' s',width=0.4)
	mlab.savefig('isosurface/P02_teor'+'%.3i'%i+'.png')

mlab.show()
