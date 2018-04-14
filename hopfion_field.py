import sympy as sym
import numpy as np
import h5py as h5
import time

p = 2
q = 3

x, y, z, t = sym.symbols('x y z t')

a = (x*x+y*y+z*z-t*t-1+2*1j*z)/(x*x+y*y+z*z-(t-1j)*(t-1j))
b = 2*(x-1j*y)/(x*x+y*y+z*z-(t-1j)*(t-1j))

Da = np.array([sym.diff(a,x),sym.diff(a,y),sym.diff(a,z)])
Db = np.array([sym.diff(b,x),sym.diff(b,y),sym.diff(b,z)])

F = p*q * a**(p-1)*b**(q-1) * np.cross(Da,Db)

N = 100
D = 1./10

#E = np.array([[[ (complex(F[0].subs(x,(i-0.5*N)/scale).subs(y,(j-0.5*N)/scale).subs(z,(k-0.5*N)/scale).subs(t,-1.5)).real)**2 + (complex(F[1].subs(x,(i-0.5*N)/scale).subs(y,(j-0.5*N)/scale).subs(z,(k-0.5*N)/scale).subs(t,-1.5)).real)**2 + (complex(F[2].subs(x,(i-0.5*N)/scale).subs(y,(j-0.5*N)/scale).subs(z,(k-0.5*N)/scale).subs(t,-1.5)).real)**2 for k in range(N+1)] for j in range(N+1)] for i in range(N+1)])

t0 = time.time()

Ex = np.array([[[complex(F[0].subs(x,D*(i-0.5*N)).subs(y,D*(j+0.5-0.5*N)).subs(z,D*(k+0.5-0.5*N)).subs(t,-1.5-0.5*D)).real for k in range(N+1)] for j in range(N+1)] for i in range(N+1)])
Ey = np.array([[[complex(F[1].subs(x,D*(i+0.5-0.5*N)).subs(y,D*(j-0.5*N)).subs(z,D*(k+0.5-0.5*N)).subs(t,-1.5-0.5*D)).real for k in range(N+1)] for j in range(N+1)] for i in range(N+1)])
Ez = np.array([[[complex(F[2].subs(x,D*(i+0.5-0.5*N)).subs(y,D*(j+0.5-0.5*N)).subs(z,D*(k-0.5*N)).subs(t,-1.5-0.5*D)).real for k in range(N+1)] for j in range(N+1)] for i in range(N+1)])

Hx = np.array([[[complex(F[0].subs(x,D*(i+0.5-0.5*N)).subs(y,D*(j-0.5*N)).subs(z,D*(k-0.5*N)).subs(t,-1.5)).imag for k in range(N+1)] for j in range(N+1)] for i in range(N+1)])
Hy = np.array([[[complex(F[1].subs(x,D*(i-0.5*N)).subs(y,D*(j+0.5-0.5*N)).subs(z,D*(k-0.5*N)).subs(t,-1.5)).imag for k in range(N+1)] for j in range(N+1)] for i in range(N+1)])
Hz = np.array([[[complex(F[2].subs(x,D*(i-0.5*N)).subs(y,D*(j-0.5*N)).subs(z,D*(k+0.5-0.5*N)).subs(t,-1.5)).imag for k in range(N+1)] for j in range(N+1)] for i in range(N+1)])

print '%.2f' % (time.time()-t0) , ' s'

fout = h5.File('hopfion_init.h5','w')
fout2 = h5.File('hopfion_init_security.h5','w')

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

fout.close()
fout2.close()
