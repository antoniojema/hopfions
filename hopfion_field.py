import numpy as np
import time
import h5py as h5
from F25 import *

N = 200
D = 0.05

t0 = time.time()
I = np.array([[[i for k in range(N+1)] for j in range(N+1)] for i in range(N+1)])
J = np.array([[[j for k in range(N+1)] for j in range(N+1)] for i in range(N+1)])
K = np.array([[[k for k in range(N+1)] for j in range(N+1)] for i in range(N+1)])

Ex = ( Fx(   D*(I-0.5*N)   , D*(J+0.5-0.5*N) , D*(K+0.5-0.5*N) , -1.5-0.5*D ) ).real
Ey = ( Fy( D*(I+0.5-0.5*N) ,   D*(J-0.5*N)   , D*(K+0.5-0.5*N) , -1.5-0.5*D ) ).real
Ez = ( Fz( D*(I+0.5-0.5*N) , D*(J+0.5-0.5*N) ,   D*(K-0.5*N)   , -1.5-0.5*D ) ).real

Hx = ( Fx( D*(I+0.5-0.5*N) ,   D*(J-0.5*N)   ,   D*(K-0.5*N)   ,    -1.5    ) ).imag
Hy = ( Fy(   D*(I-0.5*N)   , D*(J+0.5-0.5*N) ,   D*(K-0.5*N)   ,    -1.5    ) ).imag
Hz = ( Fz(   D*(I-0.5*N)   ,   D*(J-0.5*N)   , D*(K+0.5-0.5*N) ,    -1.5    ) ).imag
print time.time()-t0 , ' s'

fout = h5.File('hopfion_init.h5','w')
fout2 = h5.File('hopfion_init_sec.h5','w')

fout['Ex'] = Ex
fout['Ey'] = Ey
fout['Ez'] = Ez
fout['Hx'] = Hx
fout['Hy'] = Hy
fout['Hz'] = Hz

fout2['Ex'] = Ex
fout2['Ey'] = Ey
fout2['Ez'] = Ez
fout2['Hx'] = Hx
fout2['Hy'] = Hy
fout2['Hz'] = Hz

fout.close(); fout2.close();


#Ex_ = [[0.25 * ( Ex[i][j][N/2] + Ex[i][j-1][N/2] + Ex[i][j][N/2-1] + Ex[i][j-1][N/2-1] ) for j in np.arange(1,N+1,1)] for i in np.arange(1,N+1,1)]
#Ey_ = [[0.25 * ( Ey[i][j][N/2] + Ey[i][j][N/2-1] + Ey[i-1][j][N/2] + Ey[i-1][j][N/2-1] ) for j in np.arange(1,N+1,1)] for i in np.arange(1,N+1,1)]
#Ez_ = [[0.25 * ( Ez[i][j][N/2] + Ez[i-1][j][N/2] + Ez[i][j-1][N/2] + Ez[i-1][j-1][N/2] ) for j in np.arange(1,N+1,1)] for i in np.arange(1,N+1,1)]

#import matplotlib.pyplot as plt
#import pylab
#x = np.array([D*(i-0.5*N) for i in range(N+1)])
#y = np.array([D*(i+0.5-0.5*N) for i in range(N+1)])
#y,x = np.meshgrid(x,y)
#plt.figure(0)
#pylab.pcolor( np.abs(Ex[:][:][N/2-1]) ); pylab.colorbar();
#plt.figure(1)
#pylab.pcolor( np.abs(Ex_[:][:]) ); pylab.colorbar();
#plt.figure(2)
#pylab.pcolor( np.abs((Fx(x,y,0.05,-1.5-0.5*D)).real) ); pylab.colorbar();
#plt.show()
