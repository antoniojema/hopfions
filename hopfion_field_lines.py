import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from F23 import *

N_ = 3
N = 1000000
D = 0.0001
alpha = 0.1

x = [[0 for j in range(N)] for i in range(N_)]
y = [[0 for j in range(N)] for i in range(N_)]
z = [[0 for j in range(N)] for i in range(N_)]
if N_ == 1:
	x[0][0] = 1.
	y[0][0] = 0.
	z[0][0] = 0.
else:
	for i in range(N_):
		x[i][0] = np.cos((-alpha*(i-N_+2))/(N_-1))
		y[i][0] = np.sin((-alpha*(i-N_+2))/(N_-1))
		z[i][0] = 0.

for i in np.arange(1,N,1):
	print 100.*i/N,' %'

	for j in range(N_):
		Ex = (Fx(x[j][i-1],y[j][i-1],z[j][i-1],-1.5)).real
		Ey = (Fy(x[j][i-1],y[j][i-1],z[j][i-1],-1.5)).real
		Ez = (Fz(x[j][i-1],y[j][i-1],z[j][i-1],-1.5)).real
		E = Ex*Ex + Ey*Ey + Ez*Ez
		x[j][i] = x[j][i-1] + D * Ex/E
		y[j][i] = y[j][i-1] + D * Ey/E
		z[j][i] = z[j][i-1] + D * Ez/E


fig = plt.figure()
ax = fig.gca(projection='3d')
for i in range(N_):
	ax.plot(x[i][:],y[i][:],z[i][:],label=str((0.9*i-1.1*(i-N_+1))/(N_-1)))
plt.show()
