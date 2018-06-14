import numpy as np
import time
import h5py as h5
from F25 import *

N = 200
D = 10./N
Dt = D/2.

t0 = time.time()
I = D * (np.array([[[i for k in range(N+1)] for j in range(N+1)] for i in range(N+1)]) - 0.5*N)
J = D * (np.array([[[j for k in range(N+1)] for j in range(N+1)] for i in range(N+1)]) - 0.5*N)
K = D * (np.array([[[k for k in range(N+1)] for j in range(N+1)] for i in range(N+1)]) - 0.5*N)

Ex = ( Fx(    I    , J+0.5*D , K+0.5*D , -1.5-0.5*Dt ) ).real
Ey = ( Fy( I+0.5*D ,    J    , K+0.5*D , -1.5-0.5*Dt ) ).real
Ez = ( Fz( I+0.5*D , J+0.5*D ,    K    , -1.5-0.5*Dt ) ).real

Hx = ( Fx( I+0.5*D ,    J    ,    K    ,    -1.5     ) ).imag
Hy = ( Fy(    I    , J+0.5*D ,    K    ,    -1.5     ) ).imag
Hz = ( Fz(    I    ,    J    , K+0.5*D ,    -1.5     ) ).imag

Ex2 = ( Fx(    I    , J+0.5*D , K+0.5*D , -1.5-1.5*Dt ) ).real
Ey2 = ( Fy( I+0.5*D ,    J    , K+0.5*D , -1.5-1.5*Dt ) ).real
Ez2 = ( Fz( I+0.5*D , J+0.5*D ,    K    , -1.5-1.5*Dt ) ).real

Hx2 = ( Fx( I+0.5*D ,    J    ,    K    ,   -1.5-Dt   ) ).imag
Hy2 = ( Fy(    I    , J+0.5*D ,    K    ,   -1.5-Dt   ) ).imag
Hz2 = ( Fz(    I    ,    J    , K+0.5*D ,   -1.5-Dt   ) ).imag

print time.time()-t0 , ' s'

fout = h5.File('init_25.h5','w')

fout['Ex'] = Ex
fout['Ey'] = Ey
fout['Ez'] = Ez
fout['Hx'] = Hx
fout['Hy'] = Hy
fout['Hz'] = Hz

fout['Ex2'] = Ex2
fout['Ey2'] = Ey2
fout['Ez2'] = Ez2
fout['Hx2'] = Hx2
fout['Hy2'] = Hy2
fout['Hz2'] = Hz2

fout.close();