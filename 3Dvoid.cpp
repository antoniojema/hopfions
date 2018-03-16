#include <iostream>
#include <cmath>
#include <fstream>
//#include <H5Cpp.h>

using namespace std;
//using namespace H5;

double A, Xe, Ye, Ze;
double B, Xm, Ym, Zm;
const int N=100;
double Ex[N+1][N+1][N+1],  Ey[N+1][N+1][N+1],  Ez[N+1][N+1][N+1];
double Hx[N+1][N+1][N+1],  Hy[N+1][N+1][N+1],  Hz[N+1][N+1][N+1];
double Ex1[N+1][N+1][N+1], Ey1[N+1][N+1][N+1], Ez1[N+1][N+1][N+1];
double Hx1[N+1][N+1][N+1], Hy1[N+1][N+1][N+1], Hz1[N+1][N+1][N+1];
double Ex2[N+1][N+1][N+1], Ey2[N+1][N+1][N+1], Ez2[N+1][N+1][N+1];
double Hx2[N+1][N+1][N+1], Hy2[N+1][N+1][N+1], Hz2[N+1][N+1][N+1];

double fe(double a, double x, double y, double z){
	return A*a + Xe*x + Ye*y + Ze*z;
}
double fm(double b, double x, double y, double z){
	return B*b + Xm*x + Ym*y + Zm*z;
}

int main(){
	int n, i, j, k, iterations=100;
	double sigma, sigmam, Dx, Dy, Dz, Dt, epsilon, mu, c, L, x;
	/*ifstream fin;
	ofstream fout;*/
	
	L=10;
	sigma = sigmam = 0;
	Dx = Dy = Dz = Dt = L/N;
	epsilon = mu = 1;
	c=1;
	
	A = (1.-sigma*Dt/(2.*epsilon))/(1.+sigma*Dt/(2.*epsilon));
	Xe = (Dt/(epsilon*Dx))/(1.+sigma*Dt/(2.*epsilon));
	Ye = (Dt/(epsilon*Dy))/(1.+sigma*Dt/(2.*epsilon));
	Ze = (Dt/(epsilon*Dz))/(1.+sigma*Dt/(2.*epsilon));
	B = (1.-sigmam*Dt/(2.*mu))/(1.+sigmam*Dt/(2.*mu));
	Xm = (Dt/(mu*Dx))/(1.+sigmam*Dt/(2.*mu));
	Ym = (Dt/(mu*Dy))/(1.+sigmam*Dt/(2.*mu));
	Zm = (Dt/(mu*Dz))/(1.+sigmam*Dt/(2.*mu));
	
	//fout.open("3dimvacio.txt");
	
	//Condiciones iniciales
	for(i=0; i<=N; i++){
		for(j=0; j<=N; j++){
			for(k=0; k<=N; k++){
				Ex[i][j][k]=0; Ey[i][j][k]=0; Ez[i][j][k]=0;
				Hx[i][j][k]=0; Hy[i][j][k]=0; Hz[i][j][k]=0;
			}
		}
	}
	for(i=0; i<=N; i++){
		for(j=0; j<=N; j++){
			for(k=0; k<=N; k++){
				Ex1[i][j][k] = Ex2[i][j][k] = Ey1[i][j][k] = Ey2[i][j][k] = Ez1[i][j][k] = Ez2[i][j][k]=0;
				Hx1[i][j][k] = Hx2[i][j][k] = Hy1[i][j][k] = Hy2[i][j][k] = Hz1[i][j][k] = Hz2[i][j][k]=0;
			}
		}
	}
	
	for(n=0; n<=iterations; n++){
		//Calculo los nuevos E
		for(i=0; i<=N-1; i++){
			for(j=0; j<=N-1; j++){
				for(k=0; k<=N-1; k++){
					Ex[i][j][k] = fe(Ex[i][j][k], 0, Hz[i][j+1][k]-Hz[i][j][k], Hy[i][j][k]-Hy[i][j][k+1]);
					Ey[i][j][k] = fe(Ey[i][j][k], Hz[i][j][k]-Hz[i+1][j][k], 0, Hx[i][j][k+1]-Hx[i][j][k]);
					Ez[i][j][k] = fe(Ez[i][j][k], Hy[i+1][j][k]-Hy[i][j][k], Hx[i][j][k]-Hx[i][j+1][k], 0);
				}
				Ex[i][j][N] = fe(Ex[i][j][N], 0, Hz[i][j+1][N]-Hz[i][j][N], Hy[i][j][N]-Hy[i][j][0]);
				Ey[i][j][N] = fe(Ey[i][j][N], Hz[i][j][N]-Hz[i+1][j][N], 0, Hx[i][j][0]-Hx[i][j][N]);
				Ez[i][j][N] = fe(Ez[i][j][N], Hy[i+1][j][N]-Hy[i][j][N], Hx[i][j][N]-Hx[i][j+1][N], 0);
			}
			for(k=0; k<=N-1; k++){
				Ex[i][N][k] = fe(Ex[i][N][k], 0, Hz[i][0][k]-Hz[i][N][k], Hy[i][N][k]-Hy[i][N][k+1]);
				Ey[i][N][k] = fe(Ey[i][N][k], Hz[i][N][k]-Hz[i+1][N][k], 0, Hx[i][N][k+1]-Hx[i][N][k]);
				Ez[i][N][k] = fe(Ez[i][N][k], Hy[i+1][N][k]-Hy[i][N][k], Hx[i][N][k]-Hx[i][0][k], 0);
			}
			Ex[i][N][N] = fe(Ex[i][N][N], 0, Hz[i][0][N]-Hz[i][N][N], Hy[i][N][N]-Hy[i][N][0]);
			Ey[i][N][N] = fe(Ey[i][N][N], Hz[i][N][N]-Hz[i+1][N][N], 0, Hx[i][N][0]-Hx[i][N][N]);
			Ez[i][N][N] = fe(Ez[i][N][N], Hy[i+1][N][N]-Hy[i][N][N], Hx[i][N][N]-Hx[i][0][N], 0);
		}
		for(j=0; j<=N-1; j++){
			for(k=0; k<=N-1;k++){
				Ex[N][j][k] = fe(Ex[N][j][k], 0, Hz[N][j+1][k]-Hz[N][j][k], Hy[N][j][k]-Hy[N][j][k+1]);
				Ey[N][j][k] = fe(Ey[N][j][k], Hz[N][j][k]-Hz[0][j][k], 0, Hx[N][j][k+1]-Hx[N][j][k]);
				Ez[N][j][k] = fe(Ez[N][j][k], Hy[0][j][k]-Hy[N][j][k], Hx[N][j][k]-Hx[N][j+1][k], 0);
			}
			Ex[N][j][N] = fe(Ex[N][j][N], 0, Hz[N][j+1][N]-Hz[N][j][N], Hy[N][j][N]-Hy[N][j][0]);
			Ey[N][j][N] = fe(Ey[N][j][N], Hz[N][j][N]-Hz[0][j][N], 0, Hx[N][j][0]-Hx[N][j][N]);
			Ez[N][j][N] = fe(Ez[N][j][N], Hy[0][j][N]-Hy[N][j][N], Hx[N][j][N]-Hx[N][j+1][N], 0);
		}
		for(k=0; k<=N-1; k++){
			Ex[N][N][k] = fe(Ex[N][N][k], 0, Hz[N][0][k]-Hz[N][N][k], Hy[N][N][k]-Hy[N][N][k+1]);
			Ey[N][N][k] = fe(Ey[N][N][k], Hz[N][N][k]-Hz[0][N][k], 0, Hx[N][N][k+1]-Hx[N][N][k]);
			Ez[N][N][k] = fe(Ez[N][N][k], Hy[0][N][k]-Hy[N][N][k], Hx[N][N][k]-Hx[N][0][k], 0);
		}
		Ex[N][N][N] = fe(Ex[N][N][N], 0, Hz[N][0][N]-Hz[N][N][N], Hy[N][N][N]-Hy[N][N][0]);
		Ey[N][N][N] = fe(Ey[N][N][N], Hz[N][N][N]-Hz[0][N][N], 0, Hx[N][N][0]-Hx[N][N][N]);
		Ez[N][N][N] = fe(Ez[N][N][N], Hy[0][N][N]-Hy[N][N][N], Hx[N][N][N]-Hx[N][0][N], 0);
		
		//E -> E1 -> E2
		for(i=0; i<=N; i++){
			for(j=0; j<=N; j++){
				for(k=0; k<=N; k++){
					Ex2[i][j][k] = Ex1[i][j][k];
					Ex1[i][j][k] = Ex[i][j][k];
					Ey2[i][j][k] = Ey1[i][j][k];
					Ey1[i][j][k] = Ey[i][j][k];
					Ez2[i][j][k] = Ez1[i][j][k];
					Ez1[i][j][k] = Ez[i][j][k];
				}
			}
		}
		
		//Hago lo propio con H
		for(i=1; i<=N; i++){
			for(j=1; j<=N; j++){
				for(k=1; k<=N; k++){
					Hx[i][j][k] = fm(Hx[i][j][k], 0, Ez[i][j-1][k]-Ez[i][j][k], Ey[i][j][k]-Ey[i][j][k-1]);
					Hy[i][j][k] = fm(Hy[i][j][k], Ez[i][j][k]-Ez[i-1][j][k], 0, Ex[i][j][k-1]-Ex[i][j][k]);
					Hz[i][j][k] = fm(Hz[i][j][k], Ey[i-1][j][k]-Ey[i][j][k], Ex[i][j][k]-Ex[i][j-1][k], 0);
				}
				//Hx[i][j][0] = fm(Hx[i][j][0], 0, Ez[i][j-1][0]-Ez[i][j][0], Ey[i][j][0]-Ey[i][j][N]);
				//Hy[i][j][0] = fm(Hy[i][j][0], Ez[i][j][0]-Ez[i-1][j][0], 0, Ex[i][j][N]-Ex[i][j][0]);
				//Hz[i][j][0] = fm(Hz[i][j][0], Ey[i-1][j][0]-Ey[i][j][0], Ex[i][j][0]-Ex[i][j-1][0], 0);
				Hx[i][j][0] = -1.*Hx2[i][j][1] + (Hx1[i][j][0]+Hx1[i][j][1]) + 0.25*(Hx1[i+1][j][0]+Hx1[i+1][j][1]+Hx1[i-1][j][0]+Hx1[i-1][j][1]+Hx1[i][j+1][0]+Hx1[i][j+1][1]+Hx1[i][j-1][0]+Hx1[i][j-1][0]-4.*Hx1[i][j][0]-4.*Hx1[i][j][1]);
				Hy[i][j][0] = -1.*Hy2[i][j][1] + (Hy1[i][j][0]+Hy1[i][j][1]) + 0.25*(Hy1[i+1][j][0]+Hy1[i+1][j][1]+Hy1[i-1][j][0]+Hy1[i-1][j][1]+Hy1[i][j+1][0]+Hy1[i][j+1][1]+Hy1[i][j-1][0]+Hy1[i][j-1][0]-4.*Hy1[i][j][0]-4.*Hy1[i][j][1]);
				Hz[i][j][0] = -1.*Hz2[i][j][1] + (Hz1[i][j][0]+Hz1[i][j][1]) + 0.25*(Hz1[i+1][j][0]+Hz1[i+1][j][1]+Hz1[i-1][j][0]+Hz1[i-1][j][1]+Hz1[i][j+1][0]+Hz1[i][j+1][1]+Hz1[i][j-1][0]+Hz1[i][j-1][0]-4.*Hz1[i][j][0]-4.*Hz1[i][j][1]);
			}
			for(k=1; k<=N; k++){
				Hx[i][0][k] = fm(Hx[i][0][k], 0, Ez[i][N][k]-Ez[i][0][k], Ey[i][0][k]-Ey[i][0][k-1]);
				Hy[i][0][k] = fm(Hy[i][0][k], Ez[i][0][k]-Ez[i-1][0][k], 0, Ex[i][0][k-1]-Ex[i][0][k]);
				Hz[i][0][k] = fm(Hz[i][0][k], Ey[i-1][0][k]-Ey[i][0][k], Ex[i][0][k]-Ex[i][N][k], 0);
			}
			Hx[i][0][0] = fm(Hx[i][0][0], 0, Ez[i][N][0]-Ez[i][0][0], Ey[i][0][0]-Ey[i][0][N]);
			Hy[i][0][0] = fm(Hy[i][0][0], Ez[i][0][0]-Ez[i-1][0][0], 0, Ex[i][0][N]-Ex[i][0][0]);
			Hz[i][0][0] = fm(Hz[i][0][0], Ey[i-1][0][0]-Ey[i][0][0], Ex[i][0][0]-Ex[i][N][0], 0);
		}
		for(j=1; j<=N; j++){
			for(k=1; k<=N;k++){
				Hx[0][j][k] = fm(Hx[0][j][k], 0, Ez[0][j-1][k]-Ez[0][j][k], Ey[0][j][k]-Ey[0][j][k-1]);
				Hy[0][j][k] = fm(Hy[0][j][k], Ez[0][j][k]-Ez[N][j][k], 0, Ex[0][j][k-1]-Ex[0][j][k]);
				Hz[0][j][k] = fm(Hz[0][j][k], Ey[N][j][k]-Ey[0][j][k], Ex[0][j][k]-Ex[0][j-1][k], 0);
			}
			Hx[0][j][0] = fm(Hx[0][j][0], 0, Ez[0][j-1][0]-Ez[0][j][0], Ey[0][j][0]-Ey[0][j][N]);
			Hy[0][j][0] = fm(Hy[0][j][0], Ez[0][j][0]-Ez[N][j][0], 0, Ex[0][j][N]-Ex[0][j][0]);
			Hz[0][j][0] = fm(Hz[0][j][0], Ey[N][j][0]-Ey[0][j][0], Ex[0][j][0]-Ex[0][j-1][0], 0);
		}
		for(k=1; k<=N; k++){
			Hx[0][0][k] = fm(Hx[0][0][k], 0, Ez[0][N][k]-Ez[0][0][k], Ey[0][0][k]-Ey[0][0][k-1]);
			Hy[0][0][k] = fm(Hy[0][0][k], Ez[0][0][k]-Ez[N][0][k], 0, Ex[0][0][k-1]-Ex[0][0][k]);
			Hz[0][0][k] = fm(Hz[0][0][k], Ey[N][0][k]-Ey[0][0][k], Ex[0][0][k]-Ex[0][N][k], 0);
		}
		Hx[0][0][0] = fm(Hx[0][0][0], 0, Ez[0][N][0]-Ez[0][0][0], Ey[0][0][0]-Ey[0][0][N]);
		Hy[0][0][0] = fm(Hy[0][0][0], Ez[0][0][0]-Ez[N][0][0], 0, Ex[0][0][N]-Ex[0][0][0]);
		Hz[0][0][0] = fm(Hz[0][0][0], Ey[N][0][0]-Ey[0][0][0], Ex[0][0][0]-Ex[0][N][0], 0);
		
		//H -> H1 -> H2
		for(i=0; i<=N; i++){
			for(j=0; j<=N; j++){
				for(k=0; k<=N; k++){
					Hx2[i][j][k] = Hx1[i][j][k];
					Hx1[i][j][k] = Hx[i][j][k];
					Hy2[i][j][k] = Hy1[i][j][k];
					Hy1[i][j][k] = Hy[i][j][k];
					Hz2[i][j][k] = Hz1[i][j][k];
					Hz1[i][j][k] = Hz[i][j][k];
				}
			}
		}
		
		//Ploteo los datos en el h5
		for(i=0; i<=N; i++){
			for(j=0; j<=N; j++){
				for(k=0; k<=N; k++){
					//Es demasiado pesado, hay que hacerlo en h5 por cojones
				}
			}
		}
	}
	//fout.close();
	
	return 0;
}
