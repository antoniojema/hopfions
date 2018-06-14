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

def ind(i,j,k):
	global N
	return (N+1)*(N+1)*i + (N+1)*j + k

fin = h5.File('results_23.h5','r')

N = fin.attrs['N'][0]
iterations = fin.attrs['iterations'][0]
D = fin.attrs['Dx'][0]
Dt = fin.attrs['Dt'][0]

J, I = np.meshgrid(range(N+1),range(N+1))
I_teor, J_teor = D * ([I,J] - 0.5*N)
J_teor = D * (J - 0.5*N)
J_P, I_P = np.meshgrid( range(N) , range(N) )
I_P_teor = D * (I_P - 0.5*N)
J_P_teor = D * (J_P - 0.5*N)

### Animation ###
fig, ax = plt.subplots()

max_error = []

def init():
	global cbar, I_P, J_P, I_z, I_z, I, J, ax, field, plane, max_error
	#print 0
	
	if field == 'P':
		Ex = fin['0']['Ex'][:]
		Ey = fin['0']['Ey'][:]
		Ez = fin['0']['Ez'][:]
		if plane=='XY':
			Ex = 0.5 * ( Ex[ind(I_P,J_P,N/2)] + Ex[ind(I_P+1, J_P , N/2 )] )
			Ey = 0.5 * ( Ey[ind(I_P,J_P,N/2)] + Ey[ind( I_P ,J_P+1, N/2 )] )
			Ez = 0.5 * ( Ez[ind(I_P,J_P,N/2)] + Ez[ind( I_P , J_P ,N/2+1)] )
			Ex_teor = ( Fx( I_P_teor+0.5*D , J_P_teor+0.5*D , 0+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
			Ey_teor = ( Fy( I_P_teor+0.5*D , J_P_teor+0.5*D , 0+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
			Ez_teor = ( Fz( I_P_teor+0.5*D , J_P_teor+0.5*D , 0+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		elif plane=='XZ':
			Ex = 0.5 * ( Ex[ind(I_P,N/2,J_P)] + Ex[ind(I_P+1, N/2 , J_P )] )
			Ey = 0.5 * ( Ey[ind(I_P,N/2,J_P)] + Ey[ind( I_P ,N/2+1, J_P )] )
			Ez = 0.5 * ( Ez[ind(I_P,N/2,J_P)] + Ez[ind( I_P , N/2 ,J_P+1)] )
			Ex_teor = ( Fx( I_P_teor+0.5*D , 0+0.5*D , J_P_teor+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
			Ey_teor = ( Fy( I_P_teor+0.5*D , 0+0.5*D , J_P_teor+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
			Ez_teor = ( Fz( I_P_teor+0.5*D , 0+0.5*D , J_P_teor+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		elif plane=='YZ':
			Ex = 0.5 * ( Ex[ind(N/2,I_P,J_P)] + Ex[ind(N/2+1, I_P , J_P )] )
			Ey = 0.5 * ( Ey[ind(N/2,I_P,J_P)] + Ey[ind( N/2 ,I_P+1, J_P )] )
			Ez = 0.5 * ( Ez[ind(N/2,I_P,J_P)] + Ez[ind( N/2 , I_P ,J_P+1)] )
			Ex_teor = ( Fx( 0+0.5*D , I_P_teor+0.5*D , J_P_teor+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
			Ey_teor = ( Fy( 0+0.5*D , I_P_teor+0.5*D , J_P_teor+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
			Ez_teor = ( Fz( 0+0.5*D , I_P_teor+0.5*D , J_P_teor+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		
		P      =      Ex * Ex      +      Ey * Ey      +      Ez * Ez
		P_teor = Ex_teor * Ex_teor + Ey_teor * Ey_teor + Ez_teor * Ez_teor
		max_val = P_teor.max()
		F = abs((P-P_teor)/max_val)
	
	elif field == 'Ex':
		E = fin['0']['Ex'][:]
		if plane=='XY':
			E = E[ind(I,J,N/2)]
			E_teor = ( Fx( I_teor , J_teor+0.5*D ,  0+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		elif plane=='XZ':
			E = E[ind(I,N/2,J)]
			E_teor = ( Fx( I_teor ,  0+0.5*D , J_teor+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		elif plane=='YZ':
			E = E[ind(N/2,I,J)]
			E_teor = ( Fx(  0 , I_teor+0.5*D , J_teor+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		
		max_val = E_teor.max()
		F = abs((E-E_teor)/max_val)
	
	elif field == 'Ey':
		E = fin['0']['Ey'][:]
		if plane=='XY':
			E = E[ind(I,J,N/2)]
			E_teor = ( Fy( I_teor+0.5*D , J_teor ,  0+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		elif plane=='XZ':
			E = E[ind(I,N/2,J)]
			E_teor = ( Fy( I_teor+0.5*D ,  0 , J_teor+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		elif plane=='YZ':
			E = E[ind(N/2,I,J)]
			E_teor = ( Fy(  0+0.5*D , I_teor , J_teor+0.5*D , -1.5-0.5*Dt+0*Dt ) ).real
		
		max_val = E_teor.max()
		F = abs((E-E_teor)/max_val)
	
	elif field == 'Ez':
		E = fin['0']['Ez'][:]
		if plane=='XY':
			E = E[ind(I,J,N/2)]
			E_teor = ( Fz( I_teor+0.5*D , J_teor+0.5*D ,  0 , -1.5-0.5*Dt+0*Dt ) ).real
		elif plane=='XZ':
			E = E[ind(I,N/2,J)]
			E_teor = ( Fz( I_teor+0.5*D ,  0+0.5*D , J_teor , -1.5-0.5*Dt+0*Dt ) ).real
		elif plane=='YZ':
			E = E[ind(N/2,I,J)]
			E_teor = ( Fz(  0+0.5*D , I_teor+0.5*D , J_teor , -1.5-0.5*Dt+0*Dt ) ).real
		
		max_val = E_teor.max()
		F = abs((E-E_teor)/max_val)
	
	max_error += [F.max()]
	
	pylab.pcolor(F)
	cbar = pylab.colorbar()
        plt.clim(0,0.1)
	ax.text(1,190,'t = '+"%.3f"%(-1.5-0.5*Dt + 0*Dt),fontsize=12,color='white')
	#ax.text(1,175,'Max. teor. = '+"%.3f"%max_val,fontsize=12,color='white')
	

def iteration(n):
	global cbar, I_P, J_P, I_z, I_z, I, J, ax, field, plane, max_error
	#print n
	
	if field == 'P':
		Ex = fin[str(n)]['Ex'][:]
		Ey = fin[str(n)]['Ey'][:]
		Ez = fin[str(n)]['Ez'][:]
		if plane=='XY':
			Ex = 0.5 * ( Ex[ind(I_P,J_P,N/2)] + Ex[ind(I_P+1, J_P , N/2 )] )
			Ey = 0.5 * ( Ey[ind(I_P,J_P,N/2)] + Ey[ind( I_P ,J_P+1, N/2 )] )
			Ez = 0.5 * ( Ez[ind(I_P,J_P,N/2)] + Ez[ind( I_P , J_P ,N/2+1)] )
			Ex_teor = ( Fx( I_P_teor+0.5*D , J_P_teor+0.5*D , 0+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
			Ey_teor = ( Fy( I_P_teor+0.5*D , J_P_teor+0.5*D , 0+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
			Ez_teor = ( Fz( I_P_teor+0.5*D , J_P_teor+0.5*D , 0+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		elif plane=='XZ':
			Ex = 0.5 * ( Ex[ind(I_P,N/2,J_P)] + Ex[ind(I_P+1, N/2 , J_P )] )
			Ey = 0.5 * ( Ey[ind(I_P,N/2,J_P)] + Ey[ind( I_P ,N/2+1, J_P )] )
			Ez = 0.5 * ( Ez[ind(I_P,N/2,J_P)] + Ez[ind( I_P , N/2 ,J_P+1)] )
			Ex_teor = ( Fx( I_P_teor+0.5*D , 0+0.5*D , J_P_teor+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
			Ey_teor = ( Fy( I_P_teor+0.5*D , 0+0.5*D , J_P_teor+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
			Ez_teor = ( Fz( I_P_teor+0.5*D , 0+0.5*D , J_P_teor+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		elif plane=='YZ':
			Ex = 0.5 * ( Ex[ind(N/2,I_P,J_P)] + Ex[ind(N/2+1, I_P , J_P )] )
			Ey = 0.5 * ( Ey[ind(N/2,I_P,J_P)] + Ey[ind( N/2 ,I_P+1, J_P )] )
			Ez = 0.5 * ( Ez[ind(N/2,I_P,J_P)] + Ez[ind( N/2 , I_P ,J_P+1)] )
			Ex_teor = ( Fx( 0+0.5*D , I_P_teor+0.5*D , J_P_teor+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
			Ey_teor = ( Fy( 0+0.5*D , I_P_teor+0.5*D , J_P_teor+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
			Ez_teor = ( Fz( 0+0.5*D , I_P_teor+0.5*D , J_P_teor+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		
		P      =      Ex * Ex      +      Ey * Ey      +      Ez * Ez
		P_teor = Ex_teor * Ex_teor + Ey_teor * Ey_teor + Ez_teor * Ez_teor
		max_val = P_teor.max()
		F = abs((P-P_teor)/max_val)
	
	
	elif field == 'Ex':
		E = fin[str(n)]['Ex'][:]
		if plane=='XY':
			E = E[ind(I,J,N/2)]
			E_teor = ( Fx( I_teor , J_teor+0.5*D ,  0+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		elif plane=='XZ':
			E = E[ind(I,N/2,J)]
			E_teor = ( Fx( I_teor ,  0+0.5*D , J_teor+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		elif plane=='YZ':
			E = E[ind(N/2,I,J)]
			E_teor = ( Fx(  0 , I_teor+0.5*D , J_teor+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		
		max_val = E_teor.max()
		F = abs((E-E_teor)/max_val)
        
	elif field == 'Ey':
		E = fin[str(n)]['Ey'][:]
		if plane=='XY':
			E = E[ind(I,J,N/2)]
			E_teor = ( Fy( I_teor+0.5*D , J_teor ,  0+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		elif plane=='XZ':
			E = E[ind(I,N/2,J)]
			E_teor = ( Fy( I_teor+0.5*D ,  0 , J_teor+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		elif plane=='YZ':
			E = E[ind(N/2,I,J)]
			E_teor = ( Fy(  0+0.5*D , I_teor , J_teor+0.5*D , -1.5-0.5*Dt+n*Dt ) ).real
		
		max_val = E_teor.max()
		F = abs((E-E_teor)/max_val)
        
	elif field == 'Ez':
		E = fin[str(n)]['Ez'][:]
		if plane=='XY':
			E = E[ind(I,J,N/2)]
			E_teor = ( Fz( I_teor+0.5*D , J_teor+0.5*D ,  0 , -1.5-0.5*Dt+n*Dt ) ).real
		elif plane=='XZ':
			E = E[ind(I,N/2,J)]
			E_teor = ( Fz( I_teor+0.5*D ,  0+0.5*D , J_teor , -1.5-0.5*Dt+n*Dt ) ).real
		elif plane=='YZ':
			E = E[ind(N/2,I,J)]
			E_teor = ( Fz(  0+0.5*D , I_teor+0.5*D , J_teor , -1.5-0.5*Dt+n*Dt ) ).real
		
		max_val = E_teor.max()
		F = abs((E-E_teor)/max_val)
	
	
	max_error += [F.max()]
	
	plt.cla()
	pylab.pcolor(F)
	
	cbar.remove()
	cbar = pylab.colorbar()
        plt.clim(0,0.1)
	
	ax.text(1,190,'t = '+'%.3f'%(-1.5-0.5*Dt + n*Dt),fontsize=12,color='white')
	#ax.text(1,175,'Max. teor. = '+'%.3f'%max_val,fontsize=12,color='white')
	
	
	if n in [1,30,60,90,120,150,180,200]:
		plt.savefig('../hopfions_memoria/media/error/23_Error_'+field+'_'+plane+'_'+str(n)+'.png')


animation = ani.FuncAnimation(fig, iteration, [1,30,60,90,120,150,180,200], interval=25, init_func=init)
#plt.show()

Writer = ani.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
animation.save('error/23_Error_'+field+'_'+plane+'.mp4',writer=writer)

plt.figure()
plt.plot(max_error)
plt.xlabel('Iteration')
plt.ylabel('Max. error')
plt.savefig('error/23_Error_'+field+'_'+plane+'.png')

fout = h5.File('error/23_Error_'+field+'_'+plane+'.h5','w')
fout['max_error'] = max_error
fout.close()