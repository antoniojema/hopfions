import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.animation as ani
from mpl_toolkits.mplot3d import Axes3D
import h5py as h5
import time

def ind(i,j,k):
	global N
	return (N+1)*(N+1)*i+(N+1)*j+k

fin = h5.File('simulation_results.h5','r')

N = fin.attrs['N'][0]
#iterations = fin.attrs['iterations'][0]
Nlines = 1
iters = 250
alpha=0.1

E = np.zeros(3)
Ex = fin['0']['Ex'][:]
Ey = fin['0']['Ey'][:]
Ez = fin['0']['Ez'][:]

print "Begins"
t0 = time.time()
I = np.array( [[[i for k in range(N)] for j in range(N)] for i in range(N)] )
J = np.array( [[[j for k in range(N)] for j in range(N)] for i in range(N)] )
K = np.array( [[[k for k in range(N)] for j in range(N)] for i in range(N)] )
print time.time()-t0,' s'
t1 = time.time()
Ex_ = 0.5 * (Ex[ind(I,J,K)] + Ex[ind(I+1,J,K)])
Ey_ = 0.5 * (Ey[ind(I,J,K)] + Ey[ind(I,J+1,K)])
Ez_ = 0.5 * (Ez[ind(I,J,K)] + Ez[ind(I,J,K+1)])
del I,J,K
print time.time()-t1,' s'
print 'Total: ',time.time()-t0,'  s'

'''
print "Begins"
t0 = time.time()
Ex_ = np.array([[[0.5 * ( Ex[ind(i+1,j,k)] + Ex[ind(i,j,k)] ) for k in np.arange(1,N,1)] for j in np.arange(1,N,1)] for i in np.arange(1,N,1)])
print time.time()-t0,' s'
t1 = time.time()
Ey_ = np.array([[[0.5 * ( Ey[ind(i,j+1,k)] + Ey[ind(i,j,k)] ) for k in np.arange(1,N,1)] for j in np.arange(1,N,1)] for i in np.arange(1,N,1)])
print time.time()-t1,' s'
t2 = time.time()
Ez_ = np.array([[[0.5 * ( Ez[ind(i,j,k+1)] + Ez[ind(i,j,k)] ) for k in np.arange(1,N,1)] for j in np.arange(1,N,1)] for i in np.arange(1,N,1)])
print time.time()-t2,' s'
print 'Total: ',time.time()-t0,' s'
del t0,t1,t2
'''

#P = Ex_*Ex_+Ey_*Ey_+Ez_*Ez_
#argmax = P.argmax()
#imax = int( argmax / ((N+1)*(N+1)) )
#jmax = int( ( argmax % ((N+1)*(N+1)) ) / (N+1) )
#kmax = int( ( argmax % ((N+1)*(N+1)) ) % (N+1) )

#del P,Ex_,Ey_,Ez_

I = [[0 for j in range(iters)] for i in range(Nlines)]
J = [[0 for j in range(iters)] for i in range(Nlines)]
K = [[0 for j in range(iters)] for i in range(Nlines)]
#I[0] = imax; J[0]=jmax; K[0]=kmax;
if Nlines==1:
	I[0][0] = 3*N/5
	J[0][0] = N/2
	K[0][0] = N/2
else:
	for i in range(Nlines):
		I[i][0] = N/2 + N/10 * np.cos((-alpha*(i-Nlines+2))/(Nlines-1))
		J[i][0] = N/2 + N/10 * np.sin((-alpha*(i-Nlines+2))/(Nlines-1))
		K[i][0] = N/2

for n in np.arange(1,iters,1):
	print 100.*n/iters,' %'
	
	for i in range(Nlines):
		print i
		print n
		print I[i][n-1],' ',J[i][n-1],' ',I[i][n-1]
		print len(Ex_)
		E[0] = Ex_[I[i][n-1],J[i][n-1],K[i][n-1]]
		E[1] = Ey_[I[i][n-1],J[i][n-1],K[i][n-1]]		
		E[2] = Ez_[I[i][n-1],J[i][n-1],K[i][n-1]]
		argmax = (abs(E)).argmax()
		if E[argmax] > 0:
			if argmax==0:
				I[i][n] = I[i][n-1] + 1
				J[i][n] = J[i][n-1]
				K[i][n] = K[i][n-1]
			elif argmax ==1:
				I[i][n] = I[i][n-1]
				J[i][n] = J[i][n-1] + 1
				K[i][n] = K[i][n-1]
			else:
				I[i][n] = I[i][n-1]
				J[i][n] = J[i][n-1]
				K[i][n] = K[i][n-1] + 1
		elif E[argmax] < 0:
			if argmax==0:
				I[i][n] = I[i][n-1] - 1
				J[i][n] = J[i][n-1]
				K[i][n] = K[i][n-1]
			elif argmax ==1:
				I[i][n] = I[i][n-1]
				J[i][n] = J[i][n-1] - 1
				K[i][n] = K[i][n-1]
			else:
				I[i][n] = I[i][n-1]
				J[i][n] = J[i][n-1]
				K[i][n] = K[i][n-1] - 1
		else:
			I[i][n] = I[i][n-1]
			J[i][n] = J[i][n-1]
			K[i][n] = K[i][n-1]
			print "You shouldn't see this message"

fig = plt.figure()
ax = fig.gca(projection='3d')
if Nlines == 1:
	ax.plot(I[0][:],J[0][:],K[0][:])
else:
	for i in range(Nlines):
		ax.plot(I[i][:],J[i][:],K[i][:],label=str((0.9*i-1.1*(i-Nlines+1))/(Nlines-1)))
plt.show()

