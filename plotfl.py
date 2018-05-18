import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.animation as ani
from mpl_toolkits.mplot3d import Axes3D
import h5py as h5

def f(x):
	if x==0:
		return 1
	else:
		return 0

def Ex_interpolated(i,j,k):
	global Ex
	I = int(i); J = int(j); K = int(k);
	d = np.array([0 for n in range(8)])
	d[0] = np.sqrt(   (i-I)**2   + (j-(J+1/2))**2 + (k-(K+1/2))**2 )
	d[1] = np.sqrt(   (i-I)**2   + (j-(J+1/2))**2 + (k-(K-1/2))**2 )
	d[2] = np.sqrt(   (i-I)**2   + (j-(J-1/2))**2 + (k-(K+1/2))**2 )
	d[3] = np.sqrt(   (i-I)**2   + (j-(J-1/2))**2 + (k-(K-1/2))**2 )
	d[4] = np.sqrt( (i-(I+1))**2 + (j-(J+1/2))**2 + (k-(K+1/2))**2 )
	d[5] = np.sqrt( (i-(I+1))**2 + (j-(J+1/2))**2 + (k-(K-1/2))**2 )
	d[6] = np.sqrt( (i-(I+1))**2 + (j-(J-1/2))**2 + (k-(K+1/2))**2 )
	d[7] = np.sqrt( (i-(I+1))**2 + (j-(J-1/2))**2 + (k-(K-1/2))**2 )
	if 0 in d:
		d = [f(x) for x in d]
		return np.dot([ Ex[ind( I , J , K )],
				Ex[ind( I , J ,K-1)],
				Ex[ind( I ,J-1, K )],
				Ex[ind( I ,J-1,K-1)],
				Ex[ind(I+1, J , K )],
				Ex[ind(I+1, J ,K-1)],
				Ex[ind(I+1,J-1, K )],
				Ex[ind(I+1,J-1,K-1)] ] , d)
	else:
		d = 1./d
		D = 1./(np.sum(d))
		return D * np.dot([ Ex[ind( I , J , K )],
				    Ex[ind( I , J ,K-1)],
				    Ex[ind( I ,J-1, K )],
				    Ex[ind( I ,J-1,K-1)],
				    Ex[ind(I+1, J , K )],
				    Ex[ind(I+1, J ,K-1)],
				    Ex[ind(I+1,J-1, K )],
				    Ex[ind(I+1,J-1,K-1)] ] , d)

def Ey_interpolated(i,j,k):
	global Ey
	I = int(i); J = int(j); K = int(k);
	d = np.array([0 for n in range(8)])
	d[0] = np.sqrt( (i-(I+1/2))**2 +   (j-J)**2   + (k-(K+1/2))**2 )
	d[1] = np.sqrt( (i-(I+1/2))**2 +   (j-J)**2   + (k-(K-1/2))**2 )
	d[2] = np.sqrt( (i-(I-1/2))**2 +   (j-J)**2   + (k-(K+1/2))**2 )
	d[3] = np.sqrt( (i-(I-1/2))**2 +   (j-J)**2   + (k-(K-1/2))**2 )
	d[4] = np.sqrt( (i-(I+1/2))**2 + (j-(J+1))**2 + (k-(K+1/2))**2 )
	d[5] = np.sqrt( (i-(I+1/2))**2 + (j-(J+1))**2 + (k-(K-1/2))**2 )
	d[6] = np.sqrt( (i-(I-1/2))**2 + (j-(J+1))**2 + (k-(K+1/2))**2 )
	d[7] = np.sqrt( (i-(I-1/2))**2 + (j-(J+1))**2 + (k-(K-1/2))**2 )
	if 0 in d:
		d = [f(x) for x in d]
		return np.dot([ Ey[ind( I , J , K )],
			 	Ey[ind( I , J ,K-1)],
				Ey[ind(I-1, J , K )],
				Ey[ind(I-1, J ,K-1)],
				Ey[ind( I ,J+1, K )],
				Ey[ind( I ,J+1,K-1)],
				Ey[ind(I-1,J+1, K )],
				Ey[ind(I-1,J+1,K-1)] ] , d)
	else:
		d = 1./d
		D = 1./(np.sum(d))
		return D * np.dot([ Ey[ind( I , J , K )],
				    Ey[ind( I , J ,K-1)],
				    Ey[ind(I-1, J , K )],
				    Ey[ind(I-1, J ,K-1)],
				    Ey[ind( I ,J+1, K )],
				    Ey[ind( I ,J+1,K-1)],
				    Ey[ind(I-1,J+1, K )],
				    Ey[ind(I-1,J+1,K-1)] ] , d)

def Ez_interpolated(i,j,k):
	global Ez
	I = int(i); J = int(j); K = int(k);
	d = np.array([0 for n in range(8)])
	d[0] = np.sqrt( (i-(I+1/2))**2 + (j-(J+1/2))**2 +   (k-K)**2   )
	d[1] = np.sqrt( (i-(I-1/2))**2 + (j-(J+1/2))**2 +   (k-K)**2   )
	d[2] = np.sqrt( (i-(I+1/2))**2 + (j-(J-1/2))**2 +   (k-K)**2   )
	d[3] = np.sqrt( (i-(I-1/2))**2 + (j-(J-1/2))**2 +   (k-K)**2   )
	d[4] = np.sqrt( (i-(I+1/2))**2 + (j-(J+1/2))**2 + (k-(K+1))**2 )
	d[5] = np.sqrt( (i-(I-1/2))**2 + (j-(J+1/2))**2 + (k-(K+1))**2 )
	d[6] = np.sqrt( (i-(I+1/2))**2 + (j-(J-1/2))**2 + (k-(K+1))**2 )
	d[7] = np.sqrt( (i-(I-1/2))**2 + (j-(J-1/2))**2 + (k-(K+1))**2 )
	if 0  in d:
		d = [f(x) for x in d]
		return np.dot([ Ez[ind( I , J , K )],
				Ez[ind(I-1, J , K )],
				Ez[ind( I ,J-1, K )],
				Ez[ind(I-1,J-1, K )],
				Ez[ind( I , J ,K+1)],
				Ez[ind(I-1, J ,K+1)],
				Ez[ind( I ,J-1,K+1)],
				Ez[ind(I-1,J-1,K+1)] ] , d)
	else:
		d = 1./d
		D = 1./(np.sum(d))
		return D * np.dot([ Ez[ind( I , J , K )],
				    Ez[ind(I-1, J , K )],
				    Ez[ind( I ,J-1, K )],
				    Ez[ind(I-1,J-1, K )],
				    Ez[ind( I , J ,K+1)],
				    Ez[ind(I-1, J ,K+1)],
				    Ez[ind( I ,J-1,K+1)],
				    Ez[ind(I-1,J-1,K+1)] ] , d)

def ind(i,j,k):
	return (N+1)*(N+1)*i+(N+1)*j+k

#########################################################

fin = h5.File('/home/serron/Data/repo/hopfions/3D.h5','r')

N = fin.attrs['N'][0]
#iterations = fin.attrs['iterations'][0]
Nlines = 1
iters = 50000
D = 0.1
alpha=0.1

Ex = fin['0']['Ex'][:]
Ey = fin['0']['Ey'][:]
Ez = fin['0']['Ez'][:]

#Ex_ = np.array([[0.5 * ( Ex[ind(i+1,j,N/2)] + Ex[ind(i,j,N/2)] ) for j in np.arange(1,N,1)] for i in np.arange(1,N,1)])
#Ey_ = np.array([[0.5 * ( Ey[ind(i,j+1,N/2)] + Ey[ind(i,j,N/2)] ) for j in np.arange(1,N,1)] for i in np.arange(1,N,1)])
#Ez_ = np.array([[0.5 * ( Ez[ind(i,j,N/2+1)] + Ez[ind(i,j,N/2)] ) for j in np.arange(1,N,1)] for i in np.arange(1,N,1)])

#P = Ex_*Ex_+Ey_*Ey_+Ez_*Ez_
#argmax = P.argmax()
#imax = int( argmax / ((N+1)*(N+1)) )
#jmax = int( ( argmax % ((N+1)*(N+1)) ) / (N+1) )
#kmax = int( ( argmax % ((N+1)*(N+1)) ) % (N+1) )

#del P,Ex_,Ey_,Ez_

I = np.array( [[0 for j in range(iters)] for i in range(Nlines)] )
J = np.array( [[0 for j in range(iters)] for i in range(Nlines)] )
K = np.array( [[0 for j in range(iters)] for i in range(Nlines)] )
#I[0] = imax; J[0]=jmax; K[0]=kmax;
if Nlines==1:
	I[0][0] = 3*N/5
	J[0][0] = 3*N/5
	K[0][0] = N/2
else:
	for i in range(Nlines):
		I[i][0] = N/2 + N/10 * np.cos((-alpha*(i-Nlines+2))/(Nlines-1))
		J[i][0] = N/2 + N/10 * np.sin((-alpha*(i-Nlines+2))/(Nlines-1))
		K[i][0] = N/2

for n in np.arange(1,iters,1):
	print 100.*n/iters,' %'
	
	for i in range(Nlines):
		Ex_ = Ex_interpolated(I[i][n-1],J[i][n-1],K[i][n-1])
		Ey_ = Ey_interpolated(I[i][n-1],J[i][n-1],K[i][n-1])
		Ez_ = Ez_interpolated(I[i][n-1],J[i][n-1],K[i][n-1])
		E = Ex_*Ex_ + Ey_*Ey_ + Ez_*Ez_
		I[i][n] = I[i][n-1] + D * Ex_/E
		J[i][n] = J[i][n-1] + D * Ey_/E
		K[i][n] = K[i][n-1] + D * Ez_/E

fig = plt.figure()
ax = fig.gca(projection='3d')
if Nlines == 1:
	ax.plot(I[0][:],J[0][:],K[0][:])
else:
	for i in range(Nlines):
		ax.plot(I[i][:],J[i][:],K[i][:],label=str((0.9*i-1.1*(i-Nlines+1))/(Nlines-1)))
plt.show()

