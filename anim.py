#Este codigo abre "1D.h5", generado por "1Dh5.py" y muestra una animacion
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import h5py as h5

#Open h5 file
fin = h5.File('1D.h5','r')

N = fin.attrs['N']
iterations = fin.attrs['iterations']
L = 10
x = np.linspace(0,L,N+1)
Ex = fin['0']['Ex'][:]

### Animation ###
fig,ax = plt.subplots()
line, = ax.plot(x,Ex)

def iteration(i):
	Ex = fin[str(i)]['Ex'][:]
	
	line.set_ydata(Ex)
	return line,

animation = ani.FuncAnimation(fig, iteration, iterations, interval=25)
plt.show()
