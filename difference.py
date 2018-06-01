import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import pylab
from F23 import *
import time

fin = h5.File('simulation_results.h5','r')

N = fin.attrs['N'][0]
iterations = fin.attrs['iterations'][0]
D = fin.attrs['Dx'][0]
Dt = fin.attrs['Dt'][0]

e = 0
I = np.array([[i for j in np.arange(e,N+1-e,1)] for i in np.arange(e,N+1-e,1)])
J = np.array([[j for j in np.arange(e,N+1-e,1)] for i in np.arange(e,N+1-e,1)])

def ind(i,j,k):
	global N
	return (N+1)*(N+1)*i+(N+1)*j+k

err_Ex = [0 for i in range(iterations)]
for n in range(iterations):
	print n
	
	Ex = fin[str(n)]['Ex'][:]
	'''
	Ey = fin[str(n)]['Ey'][:]
	Ez = fin[str(n)]['Ez'][:]
	Hx = fin[str(n)]['Hx'][:]
	Hy = fin[str(n)]['Hy'][:]
	Hz = fin[str(n)]['Hz'][:]
	'''
	
	t0 = time.time()
	
	err_Ex_total = 0
	'''
	err_Ey_total = 0
	err_Ez_total = 0
	err_Hx_total = 0
	err_Hy_total = 0
	err_Hz_total = 0
	'''
	for k in np.arange(e,N+1-e,1):
		### E ###
		Field = np.array(Ex[ind(I,J,k)])
		Field_teor = np.array( ( Fx(   D*(I-0.5*N)   , D*(J+0.5-0.5*N) , D*(k+0.5-0.5*N) , -1.5+(n-0.5)*Dt ) ).real )
		err_Ex_total += ((1.*(Field - Field_teor))**2).sum() / ((N+1-2*e)*(N-1-2*e))
		
		'''
		Field = np.array(Ey[ind(I,J,k)])
		Field_teor = np.array( ( Fy( D*(I+0.5-0.5*N) ,   D*(J-0.5*N)   , D*(k+0.5-0.5*N) , -1.5+(n-0.5)*Dt ) ).real )
		err_Ey_total += (abs(1.*(Field - Field_teor)/Field_teor)).sum() / ((N+1-2*e)*(N-1-2*e))
		
		Field = np.array(Ez[ind(I,J,k)])
		Field_teor = np.array( ( Fz( D*(I+0.5-0.5*N) , D*(J+0.5-0.5*N) ,   D*(k-0.5*N)   , -1.5+(n-0.5)*Dt ) ).real )
		err_Ez_total += (abs(1.*(Field - Field_teor)/Field_teor)).sum() / ((N+1-2*e)*(N-1-2*e))
		
		### H ###
		Field = np.array(Hx[ind(I,J,k)])
		Field_teor = np.array( ( Fx( D*(I+0.5-0.5*N) ,   D*(J-0.5*N)   ,   D*(k-0.5*N)   ,    -1.5+n*Dt    ) ).imag )
		err_Hx_total += (abs(1.*(Field - Field_teor)/Field_teor)).sum() / ((N+1-2*e)*(N-1-2*e))
		
		Field = np.array(Hy[ind(I,J,k)])
		Field_teor = np.array( ( Fy(   D*(I-0.5*N)   , D*(J+0.5-0.5*N) ,   D*(k-0.5*N)   ,    -1.5+n*Dt    ) ).imag )
		err_Hy_total += (abs(1.*(Field - Field_teor)/Field_teor)).sum() / ((N+1-2*e)*(N-1-2*e))
		
		Field = np.array(Hz[ind(I,J,k)])
		Field_teor = np.array( ( Fz(   D*(I-0.5*N)   ,   D*(J-0.5*N)   , D*(k+0.5-0.5*N) ,    -1.5+n*Dt    ) ).imag )
		err_Hz_total += (abs(1.*(Field - Field_teor)/Field_teor)).sum() / ((N+1-2*e)*(N-1-2*e))
		'''
	
	err_Ex[n] = err_Ex_total / (N+1-2*e)
	'''
	err_Ey[n] = err_Ey_total / (N+1-2*e)
	err_Ez[n] = err_Ez_total / (N+1-2*e)
	err_Hx[n] = err_Hx_total / (N+1-2*e)
	err_Hy[n] = err_Hz_total / (N+1-2*e)
	err_Hz[n] = err_Hz_total / (N+1-2*e)
	'''

plt.plot(err_Ex)
plt.show()
'''
plt.plot(err_Ey)
plt.show()
plt.plot(err_Ez)
plt.show()
plt.plot(err_Hx)
plt.show()
plt.plot(err_Hy)
plt.show()
plt.plot(err_Hz)
plt.show()
'''

