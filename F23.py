def Fx(x,y,z,t):
	return 1.*(24*(((x - 1j*y)**2*(-1 - t**2 + x**2 + y**2 + (2*1j)*z + z**2)*(-4/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((16*1j)*t)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (24*t**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((16*1j)*t**3)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (4*t**4)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (8*x**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((16*1j)*t*x**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (8*t**2*x**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (4*x**4)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((8*1j)*x*y)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (16*t*x*y)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((8*1j)*t**2*x*y)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((8*1j)*x**3*y)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((8*1j)*x*y**3)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (4*y**4)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((8*1j)*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (24*t*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((24*1j)*t**2*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (8*t**3*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((8*1j)*x**2*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (8*t*x**2*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((8*1j)*y**2*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (8*t*y**2*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((8*1j)*x*y*z**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (8*y**2*z**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((8*1j)*z**3)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (8*t*z**3)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (4*z**4)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4))/(-(-1j + t)**2 + x**2 + y**2 + z**2)**3))

def Fy(x,y,z,t):
	return 1.*(24*(((x - 1j*y)**2*(-1 - t**2 + x**2 + y**2 + (2*1j)*z + z**2)*((4*1j)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (16*t)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((24*1j)*t**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (16*t**3)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((4*1j)*t**4)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((4*1j)*x**4)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (8*x*y)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((16*1j)*t*x*y)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (8*t**2*x*y)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (8*x**3*y)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((8*1j)*y**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (16*t*y**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((8*1j)*t**2*y**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (8*x*y**3)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((4*1j)*y**4)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (8*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((24*1j)*t*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (24*t**2*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((8*1j)*t**3*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (8*x**2*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((8*1j)*t*x**2*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (8*y**2*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((8*1j)*t*y**2*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((8*1j)*x**2*z**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (8*x*y*z**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (8*z**3)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((8*1j)*t*z**3)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((4*1j)*z**4)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4))/(-(-1j + t)**2 + x**2 + y**2 + z**2)**3))

def Fz(x,y,z,t):
	return 1.*(24*(((x - 1j*y)**2*(-1 - t**2 + x**2 + y**2 + (2*1j)*z + z**2)*(((-8*1j)*x)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (24*t*x)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((24*1j)*t**2*x)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (8*t**3*x)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((8*1j)*x**3)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (8*t*x**3)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (8*y)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((24*1j)*t*y)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (24*t**2*y)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((8*1j)*t**3*y)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (8*x**2*y)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((8*1j)*t*x**2*y)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((8*1j)*x*y**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (8*t*x*y**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (8*y**3)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((8*1j)*t*y**3)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (8*x*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((16*1j)*t*x*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (8*t**2*x*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (8*x**3*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((8*1j)*y*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (16*t*y*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((8*1j)*t**2*y*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((8*1j)*x**2*y*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (8*x*y**2*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((8*1j)*y**3*z)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((8*1j)*x*z**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + (8*t*x*z**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (8*y*z**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - ((8*1j)*t*y*z**2)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 - (8*x*z**3)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4 + ((8*1j)*y*z**3)/(1 + (2*1j)*t - t**2 + x**2 + y**2 + z**2)**4))/(-(-1j + t)**2 + x**2 + y**2 + z**2)**3))
