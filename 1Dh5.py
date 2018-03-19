#Este codigo genera un .h5 con un puslo gaussiano propagandose en 1D
#El archivo que abre el .h5 y plotea la animacion es "anim.py"
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import h5py as h5

def f(x):
	return np.exp(-1.*x*x)

#Grid data
N = 2000
iterations = 10001
L = 10
Dt = Dz = 1.*L/N
sigma = 0
sigmam = 0
epsilon = 1
mu = 1
c=1

#Initial conditions
x = np.linspace(0,L,N+1)
Ex = f(2*(x-4.))
Hy = f(2*(x-4.))
Ex1 = np.linspace(0,0,2)
Ex2 = np.linspace(0,0,2)
Hy1 = np.linspace(0,0,2)
Hy2 = np.linspace(0,0,2)

A=(1.-sigma*Dt/(2.*epsilon))/(1.+sigma*Dt/(2.*epsilon))
Ze=(Dt/(epsilon*Dz))/(1.+sigma*Dt/(2.*epsilon))
B=(1.-sigmam*Dt/(2.*mu))/(1.+sigmam*Dt/(2.*mu))
Zm=(Dt/(mu*Dz))/(1.+sigmam*Dt/(2.*mu))

### Simulation ###

#Inicializate h5 file
fout = h5.File('1D.h5','w')
fout.attrs['N'] = N
fout.attrs['iterations'] = iterations

#Simulation loop
for n in range(iterations):
	#Calculate Ex
	aux = np.insert(np.delete(Hy,0),N,0)
	Ex = A*Ex + Ze*(Hy-aux)
	#MUR ABC
	Ex[N] = -1.*Ex2[1] + ((c*Dt-Dz)/(c*Dt+Dz))*(Ex[N-1]+Ex2[0]) + (2.*Dz/(c*Dt+Dz))*(Ex1[0]+Ex1[1])
	
	#Ex -> Ex1 -> Ex2
	Ex2[0]=Ex1[0]
	Ex1[0]=Ex[N]
	Ex2[1]=Ex1[1]
	Ex1[1]=Ex[N-1]
	
	#Calculate Hy
	aux = np.insert(np.delete(Ex,N),0,0)
	Hy = B*Hy + Zm*(aux-Ex)
	#MUR ABC
	Hy[0] = -1.*Hy2[1] + ((c*Dt-Dz)/(c*Dt+Dz))*(Hy[1]+Hy2[0]) + (2.*Dz/(c*Dt+Dz))*(Hy1[0]+Hy1[1])
	
	#Data -> h5
	fout.create_group(str(n))
	fout[str(n)]["Ex"] = Ex
	fout[str(n)]["Hy"] = Hy

fout.close()
